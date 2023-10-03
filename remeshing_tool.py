import bpy
import copy
from mathutils import Vector
import bmesh, math
from enum import Enum

bl_info= {
    "name" : "Ball Pivot Remeshing",
    "blender": (3,00,0),
    "category" : "Object"
}


class EdgeStatus(Enum):
    ACTIVE = 0
    INNER = 1
    BOUNDARY = 2
    
    

class Cell:
    def __init__(self,id):
        self.points = []
        self.id = id
        self.avg = None
    
    def set_avg_normal(self):
        self.avg_normal = Vector((0,0,0))
        for p in self.points:
            self.avg_normal+= p.vertex.normal
        self.avg_normal.normalize()    


class Vertex_Remesh:
    def __init__(self,v):
        self.vertex = v
        self.used = False
        self.edges =[]
    
    def co(self):
        return self.vertex.co
    
    def normal(self):
        return self.vertex.normal
    
    def print(self):
        x = []
        for e in self.edges:
            x.append((e.a.vertex.index,e.b.vertex.index))
            
        return "Vertex "+ str(self.vertex.index)+" - "+str(self.used) + " - Edges: "+ str(x)

class Edge_Remesh:
    def __init__(self, a, b, op, bc):
        self.a = a
        self.b = b
        self.opposite = op
        self.status = EdgeStatus.ACTIVE
        self.center = bc
        self.prev = None
        self.next = None
        
        
    def set_link(self,e_prev,e_next):
        self.prev = e_prev
        self.next = e_next
    
    def edge_middle(self):
        return (self.a.co()+self.b.co())/2.0
    
    def __str__(self):
        return "("+ str(self.a.vertex.index)+","+ str(self.b.vertex.index)+")"
    
    def print(self):
        return "Edge ("+str(self.a.vertex.index)+","+ str(self.b.vertex.index)+") - "+ str(self.status)+ " - Edges: "+ self.prev.__str__()+self.next.__str__()
                
class Face_Remesh:
    def __init__(self,verts):
        self.verts = verts
        self.normal =(verts[0].co()-verts[1].co()).cross(verts[0].co()-verts[2].co())
        self.normal.normalize()

    
    def bmverts(self):
        return (x.vertex for x in self.verts)
    
    def print(self):
        return "Face verts: "+ str(list(x.vertex.index for x in self.verts))
    
                
class GridCell_Remesh:
    def __init__(self, points, voxel_dim, partition_dim):
        self.cells=[]
        self.points = []
        x,y,z = voxel_dim
        x /=  (partition_dim-0.001)
        y/=(partition_dim-0.001)
        z/= (partition_dim-0.001)
        x = math.ceil(x)
        y = math.ceil(y)
        z = math.ceil(z)
        if z==0:
            z=1
        if x==0:
            x=1
        if y==0:
            y=1
        for _ in range(0,int(x*y*z)):
                self.cells.append(Cell(len(self.cells)))
        
        print(len(self.cells))
        
        min = copy.copy(points[0].co())
        for p in points[1:]:
            if min.x> p.co().x:
                min.x=p.co().x
            if min.y> p.co().y:
                min.y=p.co().y
            if min.z> p.co().z:
                min.z=p.co().z
        self.min = min
        self.rho =partition_dim
        self.x = int(x)
        self.y= int(y)
        self.z = int(z)
        for p in points:
            self.points.append(p)
            weighted_coord = (p.co()-self.min)/partition_dim
            id = math.floor(weighted_coord[0])+math.floor(weighted_coord[1])*self.x+math.floor(weighted_coord[2])*self.x*self.y
            self.cells[int(id)].points.append(p)
        for cell in self.cells:
            cell.set_avg_normal()
    
    def get_cell(self,point):
        weighted_coord = (point.co()-self.min)/self.rho
        id = math.floor(weighted_coord[0])+math.floor(weighted_coord[1])*self.x+math.floor(weighted_coord[2])*self.x*self.y
        return self.cells[int(id)]

    
    def get_neighborhood_cells(self, grid, cell):
        n_cells = []
        indices = []
        for ind_x in [0,1,-1]:
            for ind_y in [0,self.y, -self.y]:
                for ind_z in [0,self.x*self.y,-self.x*self.y]:
                    ind = cell.id+ind_x+ind_y+ind_z
                    if ind>=0 and ind < len(self.cells):
                        indices.append(ind)
        
        indices.sort()
        for ind in indices:
            n_cells. append(grid.cells[ind])
        return n_cells
        
   
REGISTER_CLASS =[EdgeStatus, Edge_Remesh, Vertex_Remesh, Face_Remesh, GridCell_Remesh]


class Remesh_Pivot(bpy.types.Operator):
    bl_idname="object.remesh_pivot"
    bl_label = "Remesh pivot"
    bl_options = {'REGISTER','UNDO'}     
         

    radius: bpy.props.FloatProperty(name="Radius", default=1, min=0.01,max=1000)
           
    def execute(self, context):
            
        obj = context.object
        mesh = obj.data
        bm = bmesh.new()
        bm.from_mesh(mesh)

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        verts = []
        for v in bm.verts:
            verts.append(Vertex_Remesh(v))

        
        grid = GridCell_Remesh(verts, obj.dimensions, self.radius)


        faces_computed = []


        face, ball_center = seed_triangle(bm,grid,self.radius)


        if not face:
            print("No seed found")
            exit(0)

        faces_computed.append(face)

        e0 = Edge_Remesh(face.verts[0],face.verts[1],face.verts[2],ball_center)
        e1 = Edge_Remesh(face.verts[1],face.verts[2],face.verts[0],ball_center)
        e2 = Edge_Remesh(face.verts[2],face.verts[0],face.verts[1],ball_center)
        e0.set_link(e2,e1)
        e1.set_link(e0,e2)
        e2.set_link(e1,e0)
        face.verts[0].edges.append(e0)
        face.verts[0].edges.append(e2)
        face.verts[1].edges.append(e0)
        face.verts[1].edges.append(e1)
        face.verts[2].edges.append(e1)
        face.verts[2].edges.append(e2)
        edges = [e0,e1,e2]
        front= [e0,e1,e2]
            
            
        e_ij = get_active_edge(front)
        while(e_ij):
            candidate_point, candidate_center =ball_pivot(e_ij,grid,self.radius)
            if candidate_point and (not candidate_point.used or on_front(candidate_point,front)):
                face = Face_Remesh((e_ij.a,candidate_point,e_ij.b))
                faces_computed.append(face)
                e_ik, e_kj= join(edges, front, e_ij,candidate_point,candidate_center)
                e_ki = find_reversed_edge_on_front(e_ik)
                e_jk = find_reversed_edge_on_front(e_kj)
                if e_ki:
                    glue(e_ki,e_ik,front)
                if e_jk:
                    glue(e_jk,e_kj,front)
            else:
                e_ij.status = EdgeStatus.BOUNDARY
            e_ij = get_active_edge(front)



        mesh.clear_geometry()

        bm1 = bmesh.new()


        for v in bm.verts:
            bm1.verts.new(v.co)
        bm1.verts.ensure_lookup_table()
        bm1.verts.index_update()

        for e in edges:
            v = []
            for vert in [e.a,e.b]:
                v.append(bm1.verts[vert.vertex.index])
            if not bm1.edges.get(v):
                bm1.edges.new(v)   
        bm1.edges.ensure_lookup_table()
        bm1.edges.index_update() 
            
        for face in faces_computed:
            v = []
            for vert in face.verts:
                v.append(bm1.verts[vert.vertex.index])
            if not bm1.faces.get(v):
                bm1.faces.new(v)
        bm1.faces.ensure_lookup_table()
        bm1.faces.index_update() 
        bm1.to_mesh(mesh)
        return {'FINISHED'}
       

def seed_triangle(bm,grid, radius):
    for cell in grid.cells:
        neighbors = grid.get_neighborhood_cells(grid, cell)
        neighbor_points = []
        for _cl in neighbors:
            for _p in _cl.points:
                neighbor_points.append(_p)
        neighbor_points = grid.points
        for p in cell.points:
            neighbor_points.sort(key = lambda e: (e.co().x - p.co().x)**2 + (e.co().y - p.co().y)**2 +(e.co().z - p.co().z)**2)
            neighbor_points.remove(p)
            for p2 in neighbor_points:
                for p3 in neighbor_points:
                    if p2==p3:
                        continue
                    face = Face_Remesh((p,p2,p3))
                    if cell.avg_normal.dot(face.normal)<0:
                        continue
                    ball_center = compute_ball_center(face,radius)
                    if ball_center and ball_is_empty(ball_center,neighbor_points,radius):
                        p.used =True
                        p2.used =True
                        p3.used =True
                        return face, ball_center 
    return None,None       
        

def compute_ball_center(face,radius):
    ab = face.verts[2].co() -face.verts[0].co()
    ac = face.verts[1].co() -face.verts[0].co() 
    ab_ac_cross = ab.cross(ac)
    ac_ab_cross = ac.cross(ab)
    div = ab_ac_cross.length_squared
    try:
        to_circumcircle_center = (ab_ac_cross.cross(ab)*ac.length_squared+ ac_ab_cross.cross(ac)* ab.length_squared)/ (2.0*div)
    except:
        return None
    
    height_squared = radius**2 -to_circumcircle_center.length_squared
    if height_squared<0:
        return None
    ball_center = face.verts[0].co() + to_circumcircle_center + face.normal * math.sqrt(height_squared)
    return ball_center
    

def ball_is_empty(center,points,radius):
    for p in points:
        if (p.co()-center).length_squared < radius**2 - 0.0001:
            return False
    return True
    

def get_active_edge(front):
    while len(front)>0:
        e = front[-1]
        if e.status == EdgeStatus.ACTIVE:
            return e
        front.pop()

def clamp(num, min, max):
    return min if num < min else max if num > max else num

def ball_pivot(edge,grid,radius):
    e_middle = edge.edge_middle()
    
    to_old_center = edge.center - e_middle
    to_old_center.normalize()
    neighbours = []
    for cell in grid.get_neighborhood_cells(grid, grid.get_cell(edge.a)):
        for p in cell.points:
            neighbours.append(p)
    #debug
    neighbours = grid.points
    count = 1
    smallest_angle = math.inf
    point_smallest_angle = None
    center_of_smallest=None
    i=0
    smallest_number =0
    for p in neighbours:
        i+=1
        face = Face_Remesh((edge.a,p,edge.b))
        
        if p.normal().dot(face.normal)<0:
            continue
        c = compute_ball_center(face,radius)
        if not c:
            continue
        
        to_new_center = c-e_middle
        to_new_center.normalize()
        if to_new_center.dot(face.normal)<0:
            
            continue
        
        skip=False
        for ee in p.edges:
            other_point = ee.b if ee.a == p else ee.a
            if ee.status == EdgeStatus.INNER and (other_point == edge.a or other_point == edge.b):
                skip = True
                break
        if not skip:
            angle = math.acos(clamp(to_old_center.dot(to_new_center),-1.0,1.0))
            if to_new_center.cross(to_old_center).dot(edge.a.co()-edge.b.co())<0:
                angle += math.pi
            if angle < smallest_angle:
                smallest_angle = angle
                point_smallest_angle = p
                center_of_smallest = c
                smallest_number = 1
    if smallest_angle != math.inf:
        if ball_is_empty(center_of_smallest,neighbours,radius):
            return point_smallest_angle, center_of_smallest
    return None,None
        
def on_front(point,front):
    for e in point.edges:
        if e.status == EdgeStatus.ACTIVE:
            return True
        if e in front:
            front.remove(e)
    return False        
        
def join(edges,front,edge,point, center):
    e_ik = Edge_Remesh(edge.a,point,edge.b,center)
    e_kj = Edge_Remesh(point,edge.b,edge.a,center)
    edges.append(e_ik)
    edges.append(e_kj)
    
    e_ik.next = e_kj
    e_ik.prev = edge.prev
    edge.prev.next=e_ik
    edge.a.edges.append(e_ik)
    
    e_kj.prev = e_ik
    e_kj.next=edge.next
    edge.next.prev = e_kj
    edge.b.edges.append(e_kj)
    
    point.used = True
    point.edges.append(e_ik)
    point.edges.append(e_kj)
    
    front.append(e_ik)
    front.append(e_kj)
    edge.status = EdgeStatus.INNER
    return e_ik, e_kj
    

def find_reversed_edge_on_front(edge):
    for ee in edge.a.edges:
        if ee.a == edge.b and ee.b ==edge.a:
            return ee
    return None

def glue(edge_reversed, edge, front):
    pass

    if edge_reversed.next == edge and edge_reversed.prev == edge and edge.prev == edge_reversed and edge.next == edge_reversed:
        edge.status == EdgeStatus.INNER
        edge_reversed.status == EdgeStatus.INNER
        if edge in front:
            front.remove(edge)
        if edge_reversed in front:
            front.remove(edge_reversed)
        return
    if edge_reversed.next == edge and edge.prev == edge_reversed:
        edge_reversed.prev.next = edge.next
        edge.next.prev = edge_reversed.prev
        edge.status == EdgeStatus.INNER
        edge_reversed.status == EdgeStatus.INNER
        if edge in front:
            front.remove(edge)
        if edge_reversed in front:
            front.remove(edge_reversed)
        return
        
    if edge_reversed.prev == edge and edge.next == edge_reversed:
        edge_reversed.next.prev = edge.prev
        edge.prev.next = edge_reversed.next
        edge.status == EdgeStatus.INNER
        edge_reversed.status == EdgeStatus.INNER
        if edge in front:
            front.remove(edge)
        if edge_reversed in front:
            front.remove(edge_reversed)

        return
    
    edge_reversed.prev.next = edge.next
    edge.next.prev= edge_reversed.prev
    edge_reversed.next.prev = edge.prev
    edge.prev.next = edge_reversed.next
    edge.status == EdgeStatus.INNER
    edge_reversed.status == EdgeStatus.INNER
    if edge in front:
        front.remove(edge)
    if edge_reversed in front:
        front.remove(edge_reversed)

def menu_func(self, context):
    self.layout.operator(Remesh_Pivot.bl_idname)

def register():
    bpy.utils.register_class(Remesh_Pivot)
    bpy.types.VIEW3D_MT_object.append(menu_func)
    
def unregister():
    bpy.utils.unregister_class(Remesh_Pivot)
    bpy.types.VIEW3D_MT_object.remove(menu_func)



if __name__=="__main__":
    register()