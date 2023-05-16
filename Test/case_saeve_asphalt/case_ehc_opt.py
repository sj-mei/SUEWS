import supy as sp
import numpy as np
import multiprocessing as mp
import pandas as pd
import geatpy as ea


class suews_ESTMehc_PrmEstimation(ea.Problem):
    def __init__(self, run_ctrl_path: str, forcing_path: str, reference_path: str, num_cpu=4) -> None:
        name = 'suews_ESTMehc_PrmEstimation'          # name for ea.Problem
        M = 1                                        # the dimension of the objective function (obj.)
        maxormins = [1]                              # 1 for minimizing the obj., -1 for maximizing the obj.
        Dim = 2                                      # the dimension of the decision variables (d.v.)
        varTypes = [0] * Dim                         # 0 for continuous d.v., 1 for discrete d.v.
        lb = [0.75, 1.6 * 1e6]                      # lower bound of d.v.
        ub = [1.5, 2.4 * 1e6]                         # upper bound of d.v.
        lbin = [1, 1]                             # 0 for excluding the lower bound, 1 for including the lower bound
        ubin = [1, 1]                             # 0 for excluding the upper bound, 1 for including the upper bound
        self.df_state_init = sp.init_supy(run_ctrl_path)
        self.df_forcing = sp.util.read_forcing(forcing_path, tstep_mod=300)
        self.num_cpu = num_cpu
        self.Ts_ref = pd.read_table(reference_path, index_col=0, parse_dates=True)["Tsurf"].values
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def evalVars(self, x):
        k = np.reshape(x[:, [0]], (-1, 1))         # a.shape = [num_pop, 1]
        cp = np.reshape(x[:, [1]], (-1, 1))

        pool = mp.Pool(processes=self.num_cpu)

        num_pop, _ = k.shape
        popID_sub = np.array_split(np.arange(num_pop), self.num_cpu)

        arg_list = []
        acc = 0
        for i in range(0, self.num_cpu):
            df_state_List = [self.df_state_init.copy() for _ in range(len(popID_sub[i]))]
            for j in range(len(popID_sub[i])):
                df_state_List[j]["k_surf"] = k[j + acc, 0]
                df_state_List[j]["cp_surf"] = cp[j + acc, 0]
            acc += len(popID_sub[i])
            arg_list.append((i, df_state_List, self.df_forcing, self.Ts_ref))
        
        res = pool.starmap(assess_ehc_performance, arg_list)
        pool.close()
        pool.join()
        res = sorted(res, key=lambda x: x["pid"])
        metric = np.concatenate([r["metric"] for r in res]).reshape(-1, 1)
        
        return metric
    

def assess_ehc_performance(pid, df_state_list, df_forcing, Ts_ref):
    rmse_list = []
    for df_state in df_state_list:
        df_output, df_state_final = sp.run_supy(df_forcing, df_state, save_state=False)
        Ts_paved = df_output.loc[1, "debug"]["Tsfc_Paved"].resample("h").mean().values
        rmse = np.sum((Ts_paved[1:] - Ts_ref)**2) / np.sum((Ts_ref - np.mean(Ts_ref))**2)
        rmse_list.append(rmse)
    return {"pid": pid, "metric" : np.array(rmse_list)}


if __name__ == "__main__":
    prob = suews_ESTMehc_PrmEstimation(run_ctrl_path="./RunControl.nml", forcing_path="./Input/Saeve_asphalt_2004_data_60.txt", reference_path="./Input/Saeve_asphalt1_2004_ESTM_Ts_data_60.txt", num_cpu=4)
    algo = ea.soea_SEGA_templet(prob, ea.Population(Encoding="RI", NIND=25), MAXGEN=50, logTras=1, trappedValue=1e-4, maxTrappedCount=10)
    res = ea.optimize(algo, verbose=True, drawLog=False, drawing=0, saveFlag=False, outputMsg=True)
    obj_V = res["ObjV"][0][0]
    k_opt, cp_opt = res["Vars"][0]
    print("Objective function value: {0} with k = {1}, cp = {2}".format(obj_V, k_opt, cp_opt))