CC=g++
#CC=/usr/bin/gcc

SRC=src
INC=include
OUT=obj

# Code optimization/debugging options
CFLAGS =  -DNDEBUG -I $(INC)
DEBUG = -O4 -L/home/$(USER)/anaconda3/lib -lpython3.9 -lcrypt -lpthread -ldl  -lutil -lm -I/home/$(USER)/anaconda3/include/python3.9 -I/home/$(USER)/anaconda3/include/python3.9 -Wno-unused-result -Wsign-compare -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O3 -ffunction-sections -pipe -isystem /home/$(USER)/anaconda3/include -fdebug-prefix-map=/opt/conda/conda-bld/python-split_1661469458365/work=/usr/local/src/conda/python-3.9.13 -fdebug-prefix-map=/home/$(USER)/anaconda3=/usr/local/src/conda-prefix -I/home/$(USER)/anaconda3/include -DNDEBUG  -O3 -Wall -L/home/$(USER)/anaconda3/lib -I/$(EBROOTPYTHON)/include/python3.9 -L$(EBROOTPYTHON)/include/python3.9

tsp: $(OUT)/main.o $(OUT)/parse.o $(OUT)/graphe.o $(OUT)/node.o $(OUT)/edge.o $(OUT)/stsp.o $(OUT)/atsp.o $(OUT)/tour.o  $(OUT)/heap.o $(OUT)/hk.o  $(OUT)/bb.o $(OUT)/hk2.o  $(OUT)/binode.o  $(OUT)/biedge.o  $(OUT)/bigraphe.o  $(OUT)/hungarianMethod.o
	$(CC) $(CFLAGS) $(OUT)/main.o  $(OUT)/parse.o $(OUT)/graphe.o $(OUT)/node.o $(OUT)/edge.o $(OUT)/stsp.o  $(OUT)/atsp.o $(OUT)/tour.o $(OUT)/heap.o $(OUT)/hk.o $(OUT)/bb.o $(OUT)/hk2.o $(OUT)/binode.o  $(OUT)/biedge.o  $(OUT)/bigraphe.o $(OUT)/hungarianMethod.o -o tsp $(DEBUG)

$(OUT)/parse.o: $(SRC)/parse.cc $(INC)/parse.h $(INC)/definitions.h 
	$(CC) $(CFLAGS) -o $(OUT)/parse.o -c $(SRC)/parse.cc $(DEBUG)

$(OUT)/main.o: $(SRC)/main.cc $(INC)/parse.h $(INC)/definitions.h 
	$(CC) $(CFLAGS) -o $(OUT)/main.o -c $(SRC)/main.cc $(DEBUG)

$(OUT)/node.o: $(SRC)/node.cc  $(INC)/node.h $(SRC)/edge.cc
	$(CC) $(CFLAGS) -o $(OUT)/node.o -c $(SRC)/node.cc $(DEBUG)

$(OUT)/graphe.o: $(SRC)/graphe.cc $(INC)/graphe.h  $(INC)/node.h $(SRC)/node.cc $(SRC)/edge.cc 
	$(CC) $(CFLAGS) -o $(OUT)/graphe.o -c $(SRC)/graphe.cc $(DEBUG)

$(OUT)/edge.o: $(SRC)/edge.cc  $(INC)/node.h $(SRC)/node.cc 
	$(CC) $(CFLAGS) -o $(OUT)/edge.o -c $(SRC)/edge.cc $(DEBUG)

$(OUT)/stsp.o: $(SRC)/stsp.cc $(INC)/stsp.h  $(INC)/graphe.h $(SRC)/graphe.cc
	$(CC) $(CFLAGS) -o $(OUT)/stsp.o -c $(SRC)/stsp.cc $(DEBUG)

$(OUT)/atsp.o: $(SRC)/atsp.cc $(INC)/atsp.h  $(INC)/graphe.h $(SRC)/graphe.cc
	$(CC) $(CFLAGS) -o $(OUT)/atsp.o -c $(SRC)/atsp.cc $(DEBUG)


$(OUT)/tour.o: $(SRC)/tour.cc $(INC)/tour.h
	$(CC) $(CFLAGS) -o $(OUT)/tour.o -c $(SRC)/tour.cc $(DEBUG)
	
	
$(OUT)/heap.o: $(SRC)/heap.cc $(INC)/heap.h
	$(CC) $(CFLAGS) -o $(OUT)/heap.o -c $(SRC)/heap.cc $(DEBUG)
	
$(OUT)/hk2.o: $(SRC)/hk2.cc $(INC)/hk2.h  $(SRC)/hk.cc $(INC)/hk.h $(SRC)/stsp.cc $(INC)/stsp.h  $(INC)/graphe.h $(SRC)/graphe.cc $(SRC)/heap.cc $(INC)/heap.h
	$(CC) $(CFLAGS) -o $(OUT)/hk2.o -c $(SRC)/hk2.cc  $(DEBUG)
	
$(OUT)/hk.o: $(SRC)/hk.cc $(INC)/hk.h  $(SRC)/hk2.cc $(INC)/hk2.h $(SRC)/stsp.cc $(INC)/stsp.h  $(INC)/graphe.h $(SRC)/graphe.cc $(SRC)/heap.cc $(INC)/heap.h
	$(CC) $(CFLAGS) -o $(OUT)/hk.o -c $(SRC)/hk.cc  $(DEBUG)
	
$(OUT)/bb.o: $(SRC)/bb.cc $(INC)/bb.h $(SRC)/hk.cc $(INC)/hk.h $(SRC)/stsp.cc $(INC)/stsp.h  $(INC)/graphe.h $(SRC)/graphe.cc $(SRC)/heap.cc $(INC)/heap.h
	$(CC) $(CFLAGS) -o $(OUT)/bb.o -c $(SRC)/bb.cc $(DEBUG)
	



$(OUT)/binode.o: $(SRC)/binode.cc  $(INC)/binode.h $(SRC)/biedge.cc
	$(CC) $(CFLAGS) -o $(OUT)/binode.o -c $(SRC)/binode.cc $(DEBUG)

$(OUT)/bigraphe.o: $(SRC)/bigraphe.cc $(INC)/bigraphe.h  $(INC)/binode.h $(SRC)/binode.cc $(SRC)/biedge.cc $(SRC)/graphe.cc $(INC)/graphe.h 
	$(CC) $(CFLAGS) -o $(OUT)/bigraphe.o -c $(SRC)/bigraphe.cc $(DEBUG)

$(OUT)/biedge.o: $(SRC)/biedge.cc  $(INC)/binode.h $(SRC)/binode.cc 
	$(CC) $(CFLAGS) -o $(OUT)/biedge.o -c $(SRC)/biedge.cc $(DEBUG)

$(OUT)/hungarianMethod.o: $(SRC)/hungarianMethod.cc  $(INC)/hungarianMethod.h $(SRC)/hungarianMethod.cc 
	$(CC) $(CFLAGS) -o $(OUT)/hungarianMethod.o -c $(SRC)/hungarianMethod.cc $(DEBUG)




clean:
	rm -f $(OUT)/*.o
	rm -f tsp
