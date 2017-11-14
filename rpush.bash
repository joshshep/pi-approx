#! /bin/bash
adir=`basename "$PWD"`
server=69.166.61.14
rsync -vra --delete --exclude=".*" --exclude=".*/" . joshua.shepherd@$server:~/$adir
