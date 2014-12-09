app: app.c
	mpicc -O2 -o app app.c profiler.c critical_path.c vector.c helpers.c

clean:
	rm critPath.out graph.png graph.dot app