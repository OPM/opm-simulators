# The opm_init_check cpp program will produce two json files which
# should describe the transmissibility graph. The structure of the
# main chunk of data is a list:
#
#
# ((i1,j1,k1) , (((i2,j2,k2), T12) , ((i3,j3,k3) , T13))),    
# ((iA,jA,kA) , (((iB,jB,kB), TAB) , ((iC,jC,kC) , TAC))),
# ....
#
#
# Cell(i1,j1,k1) is connected to cells (i2,j2,k2) and (i3,j3,k3)
# respectively, with transmissibility T12 and T13
# respectively. Furthermore cell (iA,iB,iC) is connected to cells
# (iB,jB,kB) and (iC,jC,kC) with transmissibilty TAB and TAC
# respectively.


import json




class Connection(object):
    """
    The connection class holds connection information for one cell;
    including the i,j,k of the cell and the Transmissibility of connection.
    """

    def __init__(self , i1, j1 , k1 , i2 , j2 , k2 , T):
        self.i = i2 
        self.j = j2
        self.k = k2
        self.T = T

        dx = abs(i1 - i2)
        dy = abs(j1 - j2)
        dz = abs(k1 - k2)

        if dx == 1 and dy == 0 and dz == 0:
            self.dir = "X"
        elif dx == 0 and dy == 1 and dz == 0:
            self.dir = "Y"
        elif dx == 0 and dy == 0 and dz == 1:
            self.dir = "Z"
        else:
            self.dir = "NNC"


    def __str__(self):
        return "<%d,%d,%d>(T:%g)" % (self.i , self.j , self.k , self.T)



class CellConnections(object):

    def __init__(self , i,j,k):
        self.i = i
        self.j = j
        self.k = k
        self.connection_list = []
        self.connection_map = {}
        

    def __getitem__(self , index):
        if isinstance(index,int):
            return self.connection_list[index]
        else:
            return self.connection_map[index]

    def __len__(self):
        return len(self.connection_list)

    def has_key(self , dir_key):
        return self.connection_map.has_key(dir_key)



    def addConnections(self , connections):
        for ijk,T in connections:
            new_conn = Connection( self.i , self.j , self.k , ijk[0] , ijk[1] , ijk[2] , T)
            self.connection_list.append( new_conn )
            if new_conn.dir == "NNC":
                if not self.connection_map.has_key("NNC"):
                    self.connection_map["NNC"] = []
                self.connection_map["NNC"].append( new_conn )
            else:
                self.connection_map[new_conn.dir] = new_conn


    
    def __nonzero__(self):
        if len(self.connection_list) > 0:
            return True
        else:
            return False


    def connectsWith(self, i , j , k):
        for conn in self.connection_list:
            if conn.i == i and conn.j == j and conn.k == k:
                return True
            else:
                return False


    def __str__(self):
        return "<%d,%d,%d> : %s" % (self.i , self.j , self.k , [ conn.__str__() for conn in self.connection_list ])



class TransGraph(object):

    def __init__(self , nx , ny , nz):
        self.nx = nx
        self.ny = ny
        self.nz = nz
    
        self.cell_connections = [ None ] * nx * ny * nz


    def __getitem__(self, index):
        if isinstance(index , tuple):
            g = index[0] + index[1] * self.nx + index[2]*self.nx*self.ny
        else:
            g = index

        connections = self.cell_connections[g]
        if connections is None:
            k = g / (self.nx * self.ny)
            j = (g - k * self.nx * self.ny) / self.nx
            i = g - k * self.nx * self.ny - j * self.nx
            self.cell_connections[g] = CellConnections( i,j,k )
        
        return self.cell_connections[g]
        



    def addCell(self , ijk , new_connections):
        g = ijk[0] + ijk[1] * self.nx + ijk[2]*self.nx*self.ny
        connections = self[g]
        connections.addConnections( new_connections )



    def getZTrace(self , i , j):
        trace = []
        for k in range(self.nz):
            cell = self[i,j,k]
            if cell.has_key("Z"):
                trace.append( cell["Z"].T )
            else:
                trace.append( 0 )
                
        return trace


    def getXTrace(self , j , k):
        trace = []
        for i in range(self.nx):
            cell = self[i,j,k]
            if cell.has_key("X"):
                trace.append( cell["X"].T )
            else:
                trace.append( 0 )
                
        return trace


    def getYTrace(self , i , k):
        trace = []
        for j in range(self.ny):
            cell = self[i,j,k]
            if cell.has_key("Y"):
                trace.append( cell["Y"].T )
            else:
                trace.append( 0 )

        return trace


        
    
    @staticmethod
    def load(json_file):
        with open(json_file) as fileH:
            data = json.load( fileH )
        
        dims = data["dims"]
        graph = data["graph"]
        
        trans_graph = TransGraph( dims[0] , dims[1] , dims[2])

        for cell,connections in graph:
            trans_graph.addCell( cell , connections )
        
        return trans_graph

