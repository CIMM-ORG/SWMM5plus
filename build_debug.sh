#!/bin/sh

# Update settings
TOKEN='"DebugAPI" : '
NEW_VAL="${TOKEN}$DEBUG_API"
sed -i "s#$TOKEN.*#$NEW_VAL,#" $INIT_DIR/settings.json