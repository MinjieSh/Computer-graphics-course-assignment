CC = g++
CPPFLAGS = -g -Wall -O3 -I eigen-eigen-b3f3d4950030/
assignment = raytracer

$(assignment): Model.o Light.o Sphere.o Camera.o raytracer.o
	$(CC) $(CPPFLAGS)  Model.cc Light.cc Sphere.cc Camera.cc raytracer.cc -o raytracer
	
Model.o: Model.cc Model.h
	$(CC) $(CPPFLAGS) -c Model.cc

Light.o: Light.cc Light.h
	$(CC) $(CPPFLAGS) -c Light.cc

Sphere.o: Sphere.cc Sphere.h
	$(CC) $(CPPFLAGS) -c Sphere.cc
	
Camera.o: Camera.cc Camera.h
	$(CC) $(CPPFLAGS) -c Camera.cc
	
raytracer.o: raytracer.cc
	$(CC) $(CPPFLAGS) -c raytracer.cc	
	
clean:
	rm -f $(assignment) *.o *.gch
