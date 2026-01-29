# Local Setup for daily_fit.R

This guide helps you run the `daily_fit.R` script locally with secure token storage and automated scheduling.

## Quick Start

1. **Test the script (no setup needed with Yahoo Finance):**
   ```bash
   chmod +x run_docker_local.sh
   ./run_docker_local.sh
   ```

2. **Optional: Set up environment file for customization:**
   ```bash
   cp env.example .env
   # Customize COMPUTE_SOURCE or add RIINGO_TOKEN if needed
   ```

3. **Set up automated scheduling:**
   ```bash
   ./setup_cron.sh
   ```

## Understanding Compute Source

This setup now tracks where each job execution runs. Discord messages will show:
- **Source: local** - Ran via `run_local.sh` on this machine
- **Source: docker-local** - Ran via `run_docker_local.sh` on this machine
- **Source: github** - Ran via GitHub Actions in the cloud
- **Source: unknown** - COMPUTE_SOURCE environment variable not set (troubleshooting needed)

Both local and GitHub Actions jobs run simultaneously at 3:40 PM & 3:50 PM EST. The local execution provides exact timing, while GitHub Actions serves as a backup if your machine is offline.

## Secure Token Storage (Optional)

### Data Source: Yahoo Finance vs Tiingo

**Default (No Token Required):**
- The script uses Yahoo Finance by default (`get = "stock.prices"`)
- No API token needed
- Works out of the box

**Optional Tiingo Endpoint:**
- To use Tiingo instead, you need a RIINGO_TOKEN
- Uncomment `get = "tiingo"` in `scripts/daily_fit.R` (line 42)
- Get a free token at [tiingo.com](https://www.tiingo.com/)

### Setting Up .env File (Optional)

The `.env` file is already added to `.gitignore` to keep your tokens secure.

**Setup:**
```bash
# Copy the template
cp env.example .env

# Edit with your settings
nano .env
```

**Your .env file can include:**
```bash
# Optional: Only if using Tiingo endpoint
RIINGO_TOKEN=your_actual_tiingo_token_here

# Optional: Customize compute source label
COMPUTE_SOURCE=local
```

## Cron Job Setup

### Automated Setup
```bash
./setup_cron.sh
```

This will:
- Check for your `.env` file
- Create a wrapper script for cron execution
- Provide you with cron job commands to add

### Manual Cron Setup

1. **Edit your crontab:**
   ```bash
   crontab -e
   ```

2. **Add one of these schedules (see Recommended Schedule below)**

### Recommended Schedule

The recommended schedule matches GitHub Actions timing to run just before market close:

```bash
# Run at 3:40 PM and 3:50 PM EST (Monday-Friday)
40 15 * * 1-5 /Users/hyunsungkim/git/homebuying/run_cron.sh
50 15 * * 1-5 /Users/hyunsungkim/git/homebuying/run_cron.sh
```

**Why this timing?**
- Market closes at 4:00 PM EST
- Running at 3:40 PM and 3:50 PM captures near-final market data
- Matches GitHub Actions schedule (20:40 & 20:50 UTC = 3:40 PM & 3:50 PM EST)
- You'll receive messages from both sources, showing which completed first

**Note on Daylight Saving Time:**
- These times are fixed local times (15:40 and 15:50 in 24-hour format)
- During EDT (summer), GitHub Actions runs at 4:40 PM EDT (still 20:40 UTC)
- Keep local timing consistent year-round for predictability

### Cron Schedule Format
```
* * * * * command
│ │ │ │ │
│ │ │ │ └─── Day of week (0-7, Sunday = 0 or 7)
│ │ │ └───── Month (1-12)
│ │ └─────── Day of month (1-31)
│ └───────── Hour (0-23)
└─────────── Minute (0-59)
```

## Monitoring

### View Logs
```bash
# View cron execution logs
tail -f cron.log

# View system cron logs
sudo journalctl -f | grep daily_fit
```

### Test Cron Job
```bash
# Test the wrapper script manually
./run_cron.sh

# Check if cron is running
ps aux | grep cron
```

## Troubleshooting

### Common Issues

1. **"Note: RIINGO_TOKEN not set"**
   - This is informational only - the script works fine with Yahoo Finance (default)
   - Only set RIINGO_TOKEN if you want to use Tiingo instead

2. **"Docker not found"** (when using run_docker_local.sh)
   - Install Docker: [docker.com](https://www.docker.com/)
   - Or use `run_local.sh` instead (requires R installed)

3. **"R not found"** (when using run_local.sh)
   - Install Docker and use `run_docker_local.sh` instead (recommended)
   - Or install R: `brew install r` (macOS)

4. **"Permission denied"**
   - Make scripts executable: `chmod +x *.sh`

5. **Cron not running**
   - Check cron service: `sudo systemctl status cron` (Linux)
   - On macOS, cron should run automatically

### Debug Mode
```bash
# Run with verbose output
bash -x ./run_local.sh
```

## Security Notes

- ✅ `.env` file is in `.gitignore`
- ✅ Tokens are not logged in cron output
- ✅ Scripts use relative paths for portability
- ⚠️  Make sure your `.env` file has restricted permissions: `chmod 600 .env`

## File Structure

```
homebuying/
├── .env                    # Your tokens (gitignored)
├── env.example            # Template for .env
├── run_local.sh           # Direct R execution
├── run_docker_local.sh    # Docker execution
├── setup_cron.sh          # Cron setup helper
├── run_cron.sh            # Cron wrapper (auto-generated)
├── cron.log               # Execution logs
└── scripts/daily_fit.R    # Main R script
```

## Verifying Compute Source

After setup, verify that compute source tracking works:

1. **Test local execution:**
   ```bash
   ./run_local.sh
   ```
   Check Discord for message showing "Source: local"

2. **Test Docker execution:**
   ```bash
   ./run_docker_local.sh
   ```
   Check Discord for message showing "Source: docker-local"

3. **Test GitHub Actions:**
   - Wait for scheduled run, or trigger manually via GitHub UI
   - Check Discord for message showing "Source: github"

4. **Verify cron execution:**
   - Check `cron.log` for successful runs
   - Verify Discord messages appear at scheduled times

**Expected Discord Message Format:**
```
📊 Daily Leverage Update - 2026-01-20
Return today: 0.45%
Optimal leverage for tomorrow: 2.0x
Expected annualized return: 12.3%

🖥️ Source: local
⏰ Time: 2026-01-20 15:40:23 EST
```
