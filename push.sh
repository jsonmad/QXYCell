#!/bin/bash
# push.sh — commit and push all changes to GitHub
# Usage:
#   ./push.sh                        # uses a default timestamped message
#   ./push.sh "your commit message"  # uses your message

MSG="${1:-"update $(date '+%Y-%m-%d %H:%M')"}"

git add -A
git commit -m "$MSG"
git push
