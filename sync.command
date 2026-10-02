#!/bin/bash
cd /Users/vojtechzika/Research/Sequences_2026

echo "Pulling latest changes..."
if ! git pull; then
  echo ""
  echo "Pull failed (likely a conflict). Stopping here — fix it manually before syncing."
  read -n 1
  exit 1
fi

if [ -z "$(git status --porcelain)" ]; then
  echo ""
  echo "No local changes to commit."
else
  read -p "Commit message: " msg
  git add .
  git commit -m "$msg"
fi

echo ""
echo "Pushing..."
if ! git push; then
  echo ""
  echo "Push failed. Check the error above (often auth or a remote change you need to pull first)."
  read -n 1
  exit 1
fi

echo ""
echo "Done. Press any key to close."
read -n 1
