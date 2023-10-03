# Ball Pivot Reconstruction in Blender

It is an experimental project that aims to implement the ball pivot algorithm on a mesh or point cloud.
The script can be used as add-on you can install in Blender and try it on every mesh you can come up with.


## Installation
Open the Preferences, under Add-ons there is a button "Install…" in the top left corner. Select the file remeshing_tool.py.
You can enable the add on by clicking on the check box.

## Usage
Select a manifold mesh, non-manifold mesh or a point cloud, under the Object menu search for Remesh Pivot or use F3 and type the name of the command.

## Features
* Ball pivot algorithm
* Voxelization of space

## Example
There is a blend file with some example meshes, which were used for testing.

Video: https://youtu.be/JlTUJY8_STg
Notes:
* Planar meshes or point clouds can cause some problems
* If the reconstruction is not complete, try changing the readius parameter


