#!/bin/bash
set -euo pipefail

# Alert-only pre-close health check for the local HomeBuying cron runner.

PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$PATH"

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
CRON_LOG="$ROOT_DIR/cron.log"
HEALTH_LOG="$ROOT_DIR/cron-health.log"
TZ_NAME="${HOMEBUYING_HEALTH_TZ:-America/Los_Angeles}"
TOKEN_FILE="${HOMEBUYING_HEALTH_TELEGRAM_TOKEN_FILE:-}"
CHAT_ID="${HOMEBUYING_HEALTH_TELEGRAM_CHAT_ID:-}"

log() {
    printf '%s %s\n' "$(TZ="$TZ_NAME" date '+%Y-%m-%d %H:%M:%S %Z')" "$*" >> "$HEALTH_LOG"
}

send_alert() {
    local text="$1"

    log "ALERT: $text"

    if [ "${HOMEBUYING_HEALTH_DRY_RUN:-0}" = "1" ]; then
        printf '%s\n' "$text"
        return 0
    fi

    if [ -z "$TOKEN_FILE" ] || [ -z "$CHAT_ID" ]; then
        log "telegram alert not sent: HOMEBUYING_HEALTH_TELEGRAM_TOKEN_FILE or HOMEBUYING_HEALTH_TELEGRAM_CHAT_ID is unset"
        return 1
    fi

    if [ ! -r "$TOKEN_FILE" ]; then
        log "telegram alert not sent: token file is not readable"
        return 1
    fi

    HOMEBUYING_HEALTH_TEXT="$text" \
    HOMEBUYING_HEALTH_TOKEN_FILE="$TOKEN_FILE" \
    HOMEBUYING_HEALTH_CHAT_ID="$CHAT_ID" \
    python3 - <<'PY'
import os
import urllib.parse
import urllib.request

with open(os.environ["HOMEBUYING_HEALTH_TOKEN_FILE"], "r", encoding="utf-8") as f:
    token = f.read().strip()

data = urllib.parse.urlencode({
    "chat_id": os.environ["HOMEBUYING_HEALTH_CHAT_ID"],
    "text": os.environ["HOMEBUYING_HEALTH_TEXT"],
}).encode("utf-8")

request = urllib.request.Request(
    f"https://api.telegram.org/bot{token}/sendMessage",
    data=data,
)

with urllib.request.urlopen(request, timeout=10) as response:
    response.read()
PY
}

today="$(TZ="$TZ_NAME" date '+%Y-%m-%d')"
dow="$(TZ="$TZ_NAME" date '+%u')"

if [ "$dow" -gt 5 ] && [ "${HOMEBUYING_HEALTH_FORCE:-0}" != "1" ]; then
    log "OK: weekend, market cron not expected"
    exit 0
fi

threshold_epoch="$(TZ="$TZ_NAME" date -j -f '%Y-%m-%d %H:%M:%S' "$today 12:40:00" '+%s')"
now_epoch="$(date '+%s')"

if [ "$now_epoch" -lt "$threshold_epoch" ] && [ "${HOMEBUYING_HEALTH_FORCE:-0}" != "1" ]; then
    log "OK: before HomeBuying cron window"
    exit 0
fi

problems=()

crontab_text="$(crontab -l 2>/dev/null || true)"

if ! printf '%s\n' "$crontab_text" | grep -Eq '^40[[:space:]]+12[[:space:]]+\*[[:space:]]+\*[[:space:]]+1-5[[:space:]]+/Users/lovey/git/homebuying/run_cron\.sh([[:space:]]|$)'; then
    problems+=("missing 12:40 PM HomeBuying crontab entry")
fi

if ! printf '%s\n' "$crontab_text" | grep -Eq '^50[[:space:]]+12[[:space:]]+\*[[:space:]]+\*[[:space:]]+1-5[[:space:]]+/Users/lovey/git/homebuying/run_cron\.sh([[:space:]]|$)'; then
    problems+=("missing 12:50 PM HomeBuying crontab entry")
fi

if [ ! -f "$CRON_LOG" ]; then
    problems+=("cron.log is missing")
else
    log_mtime="$(stat -f '%m' "$CRON_LOG")"
    log_mtime_label="$(TZ="$TZ_NAME" date -r "$log_mtime" '+%Y-%m-%d %H:%M:%S %Z')"

    if [ "$log_mtime" -lt "$threshold_epoch" ]; then
        problems+=("cron.log has not updated since the 12:40 PM PT run window; latest mtime is $log_mtime_label")
    fi

    latest_segment="$(tail -n 700 "$CRON_LOG" | awk '
        /Running daily_fit.R/ { segment = "" }
        { segment = segment $0 ORS }
        END { printf "%s", segment }
    ')"

    if ! printf '%s\n' "$latest_segment" | grep -Fq 'Daily fit completed!'; then
        problems+=("latest HomeBuying run segment has no completion marker")
    fi

    if printf '%s\n' "$latest_segment" | grep -Eq 'Execution halted|error while loading shared libraries|No chains finished successfully|Could not find CmdStan TBB library|^Error:'; then
        problems+=("latest HomeBuying run segment contains a fatal R/Stan/Docker marker")
    fi
fi

running_processes="$(pgrep -fl 'run_cron\.sh|run_docker_local\.sh|Rscript scripts/daily_fit\.R|docker build.*homebuying:local' || true)"

if [ -n "$running_processes" ] && [ "${HOMEBUYING_HEALTH_ALLOW_RUNNING:-0}" != "1" ]; then
    problems+=("HomeBuying process is still running near market close")
fi

if [ "${#problems[@]}" -eq 0 ]; then
    log "OK: HomeBuying cron completed before market close"
    exit 0
fi

problem_text="$(printf '%s; ' "${problems[@]}")"
problem_text="${problem_text%; }"

send_alert "HomeBuying cron health check: ALERT - $problem_text. You may not have the HomeBuying trade signal before market close."
