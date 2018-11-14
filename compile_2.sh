#/usr/bin/zsh

g++ -g -c $1
g++ -g -c $2
g++ -g -c $3
g++ -g $3 $2 $1 -o MP2D
