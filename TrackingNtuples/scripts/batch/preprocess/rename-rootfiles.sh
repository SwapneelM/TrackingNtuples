i=1
for file in /eos/home-s/smehta/trackingsim/dataset/*
    do
    	echo "Working on file $file"
        cp $file /eos/home-s/smehta/trackingsim/renamed-dataset/$i.root
        ((i++))
    done