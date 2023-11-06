gUB=$(echo "$1" | grep -oE '[0-9]+' | tail -1)
cp -f $2 ./model.pt
./tsp -filein $1 -UB $gUB
rm dgl*.pkl
rm model.pt
