CC=/usr/bin/g++-7
#CC=/usr/bin/gcc

SRC=src
INC=include
OUT=obj

# Code optimization/debugging options
CFLAGS =  -DNDEBUG -g -O -I $(INC) -fPIC
DEBUG = -O4

all: hkfactors.so

hkfactors.so: $(OUT)/hkfactors.o $(OUT)/parse.o $(OUT)/graphe.o $(OUT)/node.o $(OUT)/edge.o $(OUT)/stsp.o $(OUT)/atsp.o $(OUT)/tour.o  $(OUT)/heap.o $(OUT)/hk.o  $(OUT)/bb.o $(OUT)/hk2.o  $(OUT)/binode.o  $(OUT)/biedge.o  $(OUT)/bigraphe.o  $(OUT)/hungarianMethod.o
	$(CC) $(CFLAGS) -shared $(OUT)/hkfactors.o  $(OUT)/parse.o $(OUT)/graphe.o $(OUT)/node.o $(OUT)/edge.o $(OUT)/stsp.o  $(OUT)/atsp.o $(OUT)/tour.o $(OUT)/heap.o $(OUT)/hk.o $(OUT)/bb.o $(OUT)/hk2.o $(OUT)/binode.o  $(OUT)/biedge.o  $(OUT)/bigraphe.o $(OUT)/hungarianMethod.o -o hkfactors.so $(DEBUG)

$(OUT)/parse.o: $(SRC)/parse.cc $(INC)/parse.h $(INC)/definitions.h 
	$(CC) $(CFLAGS) -o $(OUT)/parse.o -c $(SRC)/parse.cc $(DEBUG)

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

$(OUT)/hkfactors.o: $(SRC)/hkfactors.cc $(INC)/hkfactors.h $(SRC)/hk.cc $(INC)/hk.h $(SRC)/stsp.cc $(INC)/stsp.h  $(INC)/graphe.h $(SRC)/graphe.cc $(SRC)/heap.cc
	$(CC) $(CFLAGS) -o $(OUT)/hkfactors.o -c $(SRC)/hkfactors.cc $(DEBUG)
	



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
