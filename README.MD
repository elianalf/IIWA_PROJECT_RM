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
