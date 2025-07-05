from tools.settings import *



#求最大加权独立集
def max_weight_independent_set_lp(trajectory_forest, collision_set, frame):
    if len(collision_set) is 0:
        return trajectory_forest
    
    # 创建问题实例，目标是最大化
    prob = pulp.LpProblem("Maximum_Weighted_Independent_Set", pulp.LpMaximize)
    solver = pulp.CPLEX(
        msg = False, 
        options=[
            'writemps=0',      # 不生成MPS文件
            'writeprob=0',     # 不生成其他问题文件
            'solnfile=',       # 禁止生成解决方案文件
            'log=0'            # 关闭日志文件
        ]
    )
    prob.setSolver(solver)
    pulp.LpSolverDefault.msg = 0
    
    # 创建决策变量字典，节点是否被选取
    x = pulp.LpVariable.dicts("x", [branch.branch_id for branch in trajectory_forest], cat=pulp.LpBinary)
    
   
    # 目标函数：最大化总权重
    prob += pulp.lpSum([branch.GetTrajectoryWeight_DIS_Continue(frame)[0] * x[branch.branch_id] for branch in trajectory_forest])
    
    # 添加约束：如果两个节点有冲突，则不能同时选取
    for (u, v) in collision_set:
        prob += x[u.branch_id] + x[v.branch_id] <= 1
    
    # 求解问题
    solution_status = prob.solve()

    max_weight = pulp.value(prob.objective)
    # for branch in x:
    #     print(branch, x[branch].varValue)
    filtered_forest = [branch for index, branch in enumerate(trajectory_forest) if x[branch.branch_id].varValue is not None and x[branch.branch_id].varValue > 0]

    print("[" + str(frame) + "]" + "\t" + "MWIS Solved: ", len(trajectory_forest), len(filtered_forest))

    return filtered_forest

#生成冲突集
def GenerateCollision(trajectory_forest):
    collision_set = []
    for index, branch in enumerate(trajectory_forest):           
        collision = list(filter(lambda b: len(b.targets) + len(branch.targets) > \
                                len(set([x.sid for x in b.targets] + [x.sid for x in branch.targets])) and \
                                b.branch_id is not branch.branch_id, trajectory_forest))            
        

        collision_pairs = [(branch, col) for col in collision]
        collision_set = collision_set + collision_pairs

    return collision_set


