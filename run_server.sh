#!/bin/bash

source activate smallmol_mod_site_localization_dash

# This is old for a single page
# python ./dash_main_n.py
# gunicorn -w 6 --threads=12 --worker-class=gthread -b 0.0.0.0:5000 --timeout 120 --max-requests 500 --max-requests-jitter 100 --graceful-timeout 120 dash_main_n:server --access-logfile /app/logs/access.log

#!/bin/bash
# python ./main.py
gunicorn -w 2 --threads=12 --worker-class=gthread -b 0.0.0.0:5002 --timeout 120 --max-requests 500 --max-requests-jitter 100 --graceful-timeout 120 main:app --access-logfile /app/logs/access.log