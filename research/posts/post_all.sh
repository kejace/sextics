#!/usr/bin/env bash
# Post the refereed wave-2 comments to the instance discussion threads under the gh-authenticated account.
# Usage: bash /tmp/hodge/posts/post_all.sh [lane ...]      (default: all lanes with a comment.md and no posted.json)
# Skips locked threads and threads that gained comments since the saved snapshot (prints them for review instead).
set -u
cd /home/kejace/devel/sextics || exit 1
echo "posting as: $(gh api user --jq .login)"
mkdir -p /tmp/hodge/posts/_log
lanes=("$@")
if [ ${#lanes[@]} -eq 0 ]; then
  for d in /tmp/hodge/posts/*/; do l=$(basename "$d"); [ -f "$d/comment.md" ] && [ ! -f "$d/posted.json" ] && [ "$l" != "_log" ] && lanes+=("$l"); done
fi
for lane in "${lanes[@]}"; do
  f=/tmp/hodge/posts/$lane/comment.md
  node_id=$(python3 -c 'import json,sys; print(json.load(open("discussions/index.json"))["instances"][sys.argv[1]]["nodeId"])' "$lane") || { echo "$lane: no discussion id"; continue; }
  snap=$(python3 -c 'import json,sys; d=json.load(open(sys.argv[1]))["data"]["node"]; cs=d["comments"]["nodes"]; print(len(cs)+sum(len(c["replies"]["nodes"]) for c in cs))' /tmp/hodge/discussions/$lane.json 2>/dev/null || echo -1)
  now=$(gh api graphql -f id="$node_id" -f query='query($id: ID!){ node(id:$id){ ... on Discussion { locked comments(first:50){ totalCount nodes { replies(first:50){ totalCount } } } } } }' --jq '[.data.node.comments.totalCount, ([.data.node.comments.nodes[].replies.totalCount]|add // 0), .data.node.locked] | @tsv' 2>/dev/null)
  set -- $now; total=$(( ${1:-0} + ${2:-0} )); locked=${3:-false}
  if [ "$locked" = "true" ]; then echo "$lane: LOCKED, skipping"; continue; fi
  if [ "$snap" != "-1" ] && [ "$total" != "$snap" ] && [ "${FORCE:-0}" != "1" ]; then echo "$lane: thread changed since snapshot ($snap -> $total): NOT posting (re-run with FORCE=1 after reviewing)"; continue; fi
  res=$(gh api graphql -f id="$node_id" -F body=@"$f" -f query='mutation($id: ID!, $body: String!) { addDiscussionComment(input: {discussionId: $id, body: $body}) { comment { id url author { login } } } }' 2>&1)
  echo "$res" > /tmp/hodge/posts/$lane/posted.json
  url=$(echo "$res" | python3 -c 'import json,sys; d=json.load(sys.stdin); c=d["data"]["addDiscussionComment"]["comment"]; print(c["url"], c["author"]["login"], c["id"])' 2>/dev/null)
  echo "$lane: posted -> ${url:-FAILED: $(echo "$res" | cut -c1-200)}"
  echo "$(date -u +%FT%TZ) $lane ${url:-FAILED}" >> /tmp/hodge/posts/_log/posted.log
  sleep 3
done
