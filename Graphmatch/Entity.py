class Item:
    def __init__(self,id1,id2,d) -> None:
        self.src_id = int(id1)
        self.dst_id = int(id2)
        self.d = d
    
    @property
    def key(self):
        return self.d
    
    def __str__(self):
        return "match index: "+ str(self.src_id+1) +"--"+ str(self.dst_id+1) + " distance:"+str(self.d)
    

class Edge:
    def __init__(self,start_id,end_id,dist,grad) -> None:
        self.start_id = start_id
        self.end_id = end_id
        self.dist = dist
        self.grad = grad

    def __str__(self) -> str:
        return "edge: "+str(self.start_id)+"--"+str(self.end_id)+" distance:"+str(self.dist)+" gradient:"+str(self.grad)