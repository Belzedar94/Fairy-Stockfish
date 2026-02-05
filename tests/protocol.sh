#!/bin/bash
# verify protocol implementations

error()
{
  echo "protocol testing failed on line $1"
  exit 1
}
trap 'error ${LINENO}' ERR

echo "protocol testing started"

cat << EOF > uci.exp
   spawn ./stockfish
   send "uci\\n"
   expect "default spell-chess"
   expect "uciok"
   send "quit\\n"
   expect eof
EOF

cat << EOF > xboard.exp
   spawn ./stockfish
   send "xboard\\n"
   send "protover 2\\n"
   expect "feature done=1"
   send "ping\\n"
   expect "pong"
   send "quit\\n"
   expect eof
EOF

for exp in uci.exp xboard.exp
do
  echo "Testing $exp"
  timeout 5 expect $exp > /dev/null
  rm $exp
done

echo "protocol testing OK"
