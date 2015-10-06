#!/bin/zsh

pushd $(dirname $0)
rsync -e "ssh -p 222" $HOME/projects/chipboard/* pah11816@rhskl5.uni-r.de:chipboard
popd
