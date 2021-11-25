#!/bin/sh

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1KKk8d_1nQFbM9xQjTelCmWTMfK3SA7U5' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1KKk8d_1nQFbM9xQjTelCmWTMfK3SA7U5" -O trsp_3d_set1.0000000732.tar.gz && rm -rf /tmp/cookies.txt
