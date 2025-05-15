cat src/main/resources/xifdrproject.properties | tail -n +4 | sed -r 's/\\n\\//' | sed -r 's/^(([0-9]+\.)*[0-9]*)\s*$/\n\1\n/' > Changelog.md 
