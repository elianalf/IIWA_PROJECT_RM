# Per installare e compilare :
```
$ cd <catkin_ws>/src
$ git clone https://github.com/elianalf/IIWA_PROJECT.git 
$ cd ..
$ catkin_make
```


# Come eseguire il progetto :
```
Per l'esecuzione della traiettoria:
0)Salvare i file pe_800.txt e ped_800.txt nella cartella in cui da terminale si lancia il nodo
E' consigliabile salvarli nel proprio workspace all'esterno delle cartelle /devel, /src e /build. 
Es. da ~/ws$ :

1)Lanciare il world di Gazebo and spawnare il robot
$ roslaunch iiwa_gazebo iiwa_gazebo.launch

2)Lanciare il nodo 
$ roslaunch iiwa_project clik_traj.launch 

Per l'esecuzione dell controllo d'ammettenza per i test di forza:
1)Lanciare il world di Gazebo and spawnare il robot
$ roslaunch iiwa_gazebo iiwa_gazebo.launch

2)Lanciare il nodo 
$ roslaunch iiwa_project admittance_control.launch 
```

# Note:
Nella repository sono presenti due video che mostrano l'esecuzione dei due nodi. 

adm_fz1N.mp4 mostra la risposta del robot all'applicazione di una forza lineare lungo l'asse z di circa 10 N, mentre trajectory.mp4 mostra l'esecuzione della traiettoria a f=800Hz
