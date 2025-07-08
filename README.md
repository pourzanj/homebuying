# homebuying

# Updating renv
Use `renv::snapshot(type="all")` to perform snapshots, otherwise, some vital packages such as riingo may be excluded.

## crontab
For setting up on new systems, update the riingo token in `crontab-mac` and then use the command `crontab crontab-mac`