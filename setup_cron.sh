#!/bin/bash

# Setup script for cron job
# This script helps you set up the daily_fit.R to run automatically

echo "Setting up cron job for daily_fit.R..."

# Get the current directory (absolute path)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_docker_local.sh"

echo "Script directory: $SCRIPT_DIR"
echo "Run script path: $RUN_SCRIPT"

# Check if .env file exists (optional)
if [ ! -f "$SCRIPT_DIR/.env" ]; then
    echo "Note: .env file not found in $SCRIPT_DIR"
    echo "This is OK - the script will use Yahoo Finance (no token required)."
    echo "To customize settings, you can: cp env.example .env"
    echo ""
else
    echo ".env file found - using custom settings."
fi

# Make sure the run script is executable
chmod +x "$RUN_SCRIPT"

# Create a wrapper script that sets the working directory
WRAPPER_SCRIPT="$SCRIPT_DIR/run_cron.sh"
cat > "$WRAPPER_SCRIPT" << EOF
#!/bin/bash
# Wrapper script for cron execution
# Set PATH so cron can find docker and other commands
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:\$PATH"

cd "$SCRIPT_DIR"
exec "$RUN_SCRIPT" >> "$SCRIPT_DIR/cron.log" 2>&1
EOF

chmod +x "$WRAPPER_SCRIPT"

echo ""
echo "Cron setup complete!"
echo ""
echo "Recommended cron schedules (choose one):"
echo ""
echo "NOTE: macOS cron does not support CRON_TZ. Times below are in Pacific (PT)."
echo "      Eastern equivalents shown in parentheses."
echo ""
echo "Option 1: Pre-market close at 3:40/3:50 PM ET (RECOMMENDED)"
echo "40 12 * * 1-5 $WRAPPER_SCRIPT"
echo "50 12 * * 1-5 $WRAPPER_SCRIPT"
echo ""
echo "Option 2: Post-market close at 4:05/4:15 PM ET (ensures all data available)"
echo "5 13 * * 1-5 $WRAPPER_SCRIPT"
echo "15 13 * * 1-5 $WRAPPER_SCRIPT"
echo ""
echo "Option 3: Evening run at 7:40/7:50 PM ET"
echo "40 16 * * 1-5 $WRAPPER_SCRIPT"
echo "50 16 * * 1-5 $WRAPPER_SCRIPT"
echo ""
echo "To add the cron job, run: crontab -e"
echo ""
echo "Logs will be written to: $SCRIPT_DIR/cron.log"
echo ""
echo "To test the cron job manually:"
echo "$WRAPPER_SCRIPT"
