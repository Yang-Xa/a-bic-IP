#!/usr/bin/env python
# /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *   GOBNILP Copyright (C) 2012-2017 James Cussens, Mark Bartlett        *
# *                                                                       *
# *   This program is free software; you can redistribute it and/or       *
# *   modify it under the terms of the GNU General Public License as      *
# *   published by the Free Software Foundation; either version 3 of the  *
# *   License, or (at your option) any later version.                     *
# *                                                                       *
# *   This program is distributed in the hope that it will be useful,     *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
# *   General Public License for more details.                            *
# *                                                                       *
# *   You should have received a copy of the GNU General Public License   *
# *   along with this program; if not, see                                *
# *   <http://www.gnu.org/licenses>.                                      *
"""
Simplified Python version of GOBNILP.
Only accepts precomputed scores as input in 'Jaakkola' format.
Author: James Cussens, University of York
"""
from gurobipy import *
from itertools import combinations, permutations
import sys

import argparse

parser = argparse.ArgumentParser(description='Use Gurobi for Bayesian network learning')
parser.add_argument('scores', help="Local scores in 'Jaakkola' format")
parser.add_argument('-n', '--nsols', type=int, default=1, help="How many BNs to find")
parser.add_argument('-k', '--kbest', action='store_true', help="Whether to find 'k-best' networks")
parser.add_argument('-m', '--mec', action='store_true',
                    help="Whether to only allow one DAG per Markov Equivalence Class")
# parser.add_argument('-c', '--chordal', action='store_true', help = "Whether to only allow DAGs equivalent to chordal graphs")
parser.add_argument('-v', '--verbose', action="store_true", help="Whether to be verbose")
args = parser.parse_args()


class Gobnilp(Model):
    '''Subclass of Gurobi Model class, specific 
    to learning Bayesian networks
    '''

    def __init__(self, name):
        '''Initialise a Gobnilp object
        '''
        super(Gobnilp, self).__init__(name)

        # Setting Gurobi model parameters 
        self.ModelSense = -1  # maximise objective
        self.Params.PreCrush = 1  # since (always) adding cuts
        self.Params.CutPasses = 100000  # want to allow many cuts
        self.Params.GomoryPasses = 100000  # want to allow many cuts
        self.Params.MIPFocus = 2  # focus on proving optimality
        # self.Params.Presolve = 2
        # self.Params.PreDual = 2

        # gobnilp parameters
        self._subip_cutoff = -0.99
        self._subip_cutting_timelimit = 5
        self._user_cuts_rounds_count = 0
        self._user_cuts_rounds_limit = None
        self._user_cuts_stalling_limit = None
        self._user_enforcement_rounds_count = 0
        self._user_enforcement_rounds_limit = None
        self._max_cluster_size = None

        self._cs = set()

        # remaining attributes will eventually get useful values
        # if they are needed
        self._family_scores = None
        self._bn_variables = None
        self._n = None
        self._fv_list = None
        self._child_list = None
        self._parentset_list = None
        self._get_idx = None
        self._ordered_parentsets = None

        # dictionaries for MIP variables
        self._family = None
        self._arrow = None
        self._adj = None
        self._gen = None
        self._gendiff = None
        self._absgendiff = None
        self._total_order = None

        self._c2 = None
        self._c3 = None

        self._enforcing_cluster_constraints = None
        self._enforcing_matroid_constraints = None
        self._adding_cluster_cuts = None
        self._adding_matroid_cuts = None

        self._total_order_lazy = None

        # attributes for one dag per mec constraint
        self._one_dag_per_MEC = None
        self._MEC_constraint_dynamic = None
        self._MEC_constraint_careful = None
        self._mecrep = None

        self._last_bound = None

    def dag2mec(self, dag):
        '''Return a Markov equivalence class 
        (represented by a characteristic imset - only size <= 3 
        components given ) from a DAG represented by indices
        '''
        c = set()
        # print 'MEC for'
        for i in dag:
            child = self._child_list[i]
            parent_set = self._parentset_list[i]
            # print '{0}<-{1}'.format(child, ','.join(sorted(parent_set)))
            child_singleton = frozenset([child])
            for size in range(1, len(parent_set) + 1):
                for subset in combinations(parent_set, size):
                    c.add(frozenset(subset) | child_singleton)
        # print 'is ', c
        return frozenset(c)

    def mec2dag(self, c, vs=None):
        '''Return a DAG representative from a Markov equivalence class.
        Markov equivalence class represented by a characteristic imset
        (only size <=3 components given).
        DAG represented by a list of indices.
        The representative determined by order of variables in "vs"
        (which defaults to the "natural" ordering of the BN variables).
        If the representative cannot be represented by the family variables
        in the problem then return None
        '''
        fvs = []
        if vs is None:
            vs = self._bn_variables
        vs = list(vs)  # so as not to destroy the input!
        star = {}
        for v in vs:
            star[v] = set()
        for ci in c:
            for v in ci:
                star[v].add(ci)
        while vs != []:
            for v in tuple(vs):  # work on a copy
                max_size = 0
                nbrs = set()
                biggest = frozenset()
                # find all vertices adjacent to v and
                # find biggest c-imset component containing v
                for ci in star[v]:
                    if len(ci) == 2:
                        nbrs.update(ci)
                    if len(ci) > max_size:
                        biggest = ci
                        max_size = len(ci)
                if frozenset(nbrs) == biggest:
                    # we can make all adjacent vertices into parents
                    # without alterning c-imset representation
                    fvs.append((v, biggest - frozenset([v])))
                    # remove v and the relevant c-imset components                
                    vs.remove(v)
                    for ci in star[v]:
                        for v2 in ci:
                            if v2 != v:
                                star[v2].remove(ci)
                    del star[v]
                    break
            else:
                raise ValueError("Couldn't find a sink.")
        try:
            return [self._get_idx[family] for family in fvs]
        except KeyError:
            return None

    def read_local_scores(self, fobj):
        '''Read local scores for families (variables+parent set)
        from a file.
        Assumes "Jaakkola" format. Not much error checking.
        '''
        family_scores = {}
        n = int(fobj.readline())
        if args.verbose:
            # print >>sys.stderr, 'Problem has %s variables' % n
            print('Problem has %s variables' % n, file=sys.stderr)
        fields = fobj.readline().rstrip().split()

        def init(fields):
            # no check that there are only two fields
            return fields[0], int(fields[1]), 0, {}

        current_variable, nscores, i, this_dkt = init(fields)
        for line in fobj:
            fields = line.rstrip().split()
            if i < nscores:
                # don't bother checking that fields[1] correctly specifies
                # the number of parents
                this_dkt[frozenset(fields[2:])] = float(fields[0])
                i += 1
            else:
                family_scores[current_variable] = this_dkt
                current_variable, nscores, i, this_dkt = init(fields)
        family_scores[current_variable] = this_dkt
        if args.verbose:
            # print >>sys.stderr, 'Scores read in'
            print('Scores read in', file=sys.stderr)
        ordered_parentsets = {}
        for child, scored_parentsets in family_scores.items():
            ordered_parentsets[child] = sorted(scored_parentsets, key=lambda x: scored_parentsets[x], reverse=True)
        self._family_scores = family_scores
        self._bn_variables = sorted(family_scores)
        self._n = len(self._bn_variables)
        self._ordered_parentsets = ordered_parentsets

    def add_constraints_one_dag_per_MEC(self, flag=True, dynamic=True, careful=False):
        '''With default arguments adds a constraint that only one DAG
        per Markov equivalence class is feasible.
        Constraint is effected by adding appropriate lazy constraints when Gurobi
        generates a DAG solution. If dynamic is True then which DAG is feasible for each DAG
        is arbitrary (it will be the first one Gurobi comes across). If dynamic is False
        then the representative DAG is fixed (and is determined by the method mec2dag).
        The "careful" argument is only relevant if dynamic is True. If careful is True
        then all the lazy constraints are stored (not just posted) 
        to ensure that new solutions satisfy them.
        '''
        self._one_dag_per_MEC = flag
        if flag:
            self.Params.LazyConstraints = 1
            self._MEC_constraint_dynamic = dynamic
            if dynamic:
                self._MEC_constraint_careful = careful
                if careful:
                    self._mecrep = {}

    def sol2fvs(self):
        '''Extract the family variables set to true by some Gurobi solution
        '''
        families = []
        fvs = []
        for v in self._bn_variables:
            for parent_set, var in self._family[v].items():
                if var.Xn > 0.5:
                    families.append((v, parent_set))
                    fvs.append(var)
                    break
        return families, fvs

    def enforce_MEC_representative(self):
        '''This should only be called from callback with where == GRB.Callback.MIPSOL
        i.e. when we have a solution and where it is guaranteed that 
        the solution is a DAG.
        Adds a constraint that the Markov equivalence class corresponding to the
        current solution has only one feasible DAG (its distinguished representative).

        If dynamic = True, the DAG given by the current solution will normally become the 
        distinguished representative of its Markov equivalence class.
        In this case the current solution is not cut off.
        An exception to this occurs if (1) we are being careful (self._MEC_constraint_careful is set to
        True) and (2) we have previously seen a different DAG in this MEC - but somehow
        this did not prevent the current solution being proposed. In this case the earlier
        DAG is set (again) to be the distinguished rep.

        If dynamic = False, the DAG defined by the dag2mec method is the 
        distinguished representative of the Markov equivalence class for 
        the current solution. In this case the current solution most probably will be cut off.
        '''
        dag = []
        for i, val in enumerate(self.cbGetSolution(self._fv_list)):
            if val > 0.5:  # not '== 1' due to numerical problems!
                dag.append(i)
        assert len(dag) == len(self._bn_variables)
        mec = self.dag2mec(dag)
        if self._MEC_constraint_dynamic == False:
            dag = None
        self.add_constraints_one_DAG_per_MEC(mec, dag)

    def add_constraints_one_DAG_per_MEC(self, mec, dag=None):
        '''Add the constraint that the input dag is the only feasible
        DAG in the input Markov equivalence class mec.
        If dag not supplied the distinguished representative of mec
        generated by the method mec2dag is used.
        Should only be called from within a callback where we have a solution,
        i.e. where where == GRB.Callback.MIPSOL
        The constraint is a conjunction of linear constraints.
        dag is a list of the family variables set to 1.
        If supplied there is no check that any supplied DAG is a member of the MEC mec.
        '''
        if dag is None:
            # representative DAG is fixed
            rep_dag = self.mec2dag(mec)
        elif self._MEC_constraint_careful:
            try:
                # following assignment only happens
                # if constraint for this MEC previously added
                # but somehow the supplied dag has not been ruled out
                rep_dag = self._mecrep[mec]
            except KeyError:
                # following assignment occurs if this is the first DAG
                # for this MEC which has been encountered
                rep_dag = dag
                self._mecrep[mec] = rep_dag
        else:
            # simply assume/hope that this is the first DAG seen for this MEC
            rep_dag = dag

        if self._c2 is None:
            self.add_cimset_info()
        rhs = len(self._c2) + len(self._c3) - 1
        # construct a linear expression which takes its maximum
        # value of len(self._c2) + len(self._c3) iff
        # the family variable values encode a DAG in the MEC
        # represented by mec
        lexpr = LinExpr()
        for ci, fvlist in self._c2.items() + self._c3.items():
            if ci in mec:
                lexpr.addTerms([1] * len(fvlist), fvlist)
            else:
                lexpr.addTerms([-1] * len(fvlist), fvlist)
                lexpr.addConstant(1)

        if rep_dag is None:
            # this mec not allowed since its distinguished rep
            # not available
            self.cbLazy(lexpr, GRB.LESS_EQUAL, rhs)
        else:
            for i in rep_dag:
                self.cbLazy(lexpr, GRB.LESS_EQUAL, rhs + self._fv_list[i])

    def print_simple_output(self):
        '''Print simple representation of a learned DAG to standard output
        '''
        print('**********')
        print('Optimal BN has score', self.PoolObjVal)
        print('**********')
        families, fvs = self.sol2fvs()
        for i, (child, parent_set) in enumerate(families):
            print('{0}<-{1} {2}'.format(child, ','.join(sorted(parent_set)), fvs[i].Obj))
        print('**********')

    def add_cimset_info(self):
        '''Create all c-imset components of size 2 and 3 which might possibly
        take value 1
        '''
        c2 = {}
        c3 = {}
        for child, parent_dict in self._family.items():
            for parent_set, fv in parent_dict.items():
                for parent in parent_set:
                    try:
                        c2[frozenset([child, parent])].append(fv)
                    except KeyError:
                        c2[frozenset([child, parent])] = [fv]
                for (pa1, pa2) in combinations(parent_set, 2):
                    try:
                        c3[frozenset([child, pa1, pa2])].append(fv)
                    except KeyError:
                        c3[frozenset([child, pa1, pa2])] = [fv]
        self._c2 = c2
        self._c3 = c3

    def add_variables_family(self, branch_priority=0, best_first_branch_priority=False):
        '''Adds binary "family" variables to the model
        and creates a dictionary mapping each BN variable to 
        a dictionary mapping each parent set to the relevant family variable
        The family variable corresponding to v,Pa is set to 1 iff
        there is an Pa is the set of parents of v
        All these variables have objective value set to be some 
        (previously computed) local score.
        The branching priority of these variables can also be set (default = 0)
        '''
        family = {}
        fv_list = []
        child_list = []
        parentset_list = []
        get_idx = {}
        n = 0
        for child, parent_dict in self._family_scores.items():
            family[child] = {}
            for parent_set, score in parent_dict.items():
                fv = self.addVar(
                    obj=score,
                    vtype=GRB.BINARY
                )
                fv.BranchPriority = branch_priority
                family[child][parent_set] = fv
                fv_list.append(fv)
                child_list.append(child)
                parentset_list.append(parent_set)
                get_idx[child, parent_set] = n
                n += 1
            if best_first_branch_priority:
                for i, parent_set in enumerate(self._ordered_parentsets[child]):
                    family[child][parent_set].BranchPriority = -i
        if args.verbose:
            print('%d family variables declared' % n, file=sys.stderr)
        self._family = family
        self._fv_list = fv_list
        self._child_list = child_list
        self._parentset_list = parentset_list
        self._get_idx = get_idx

    def add_variables_local_scores(self, branch_priority=0):
        '''For each variable create a continuous variable which is the local score
        of the selected parent set for that variable
        '''
        local_score = {}
        n = 0
        for child, score_dkt in self._family_scores.items():
            v = self.addVar(
                ub=max(score_dkt.values()),
                lb=min(score_dkt.values()),
                obj=1,
                vtype=GRB.CONTINUOUS)
            v.BranchPriority = branch_priority
            local_score[child] = v
            n += 1
        if args.verbose:
            print('%d local score variables declared' % n, file=sys.stderr)
        self._local_score = local_score

    def add_variables_arrow(self, branch_priority=0):
        '''Adds binary "arrow" variables to the model
        and creates a dictionary mapping each ordered pair of BN variables to its
        corresponding arrow variable.
        The arrow variable corresponding to (v1,v2) is set to 1 iff
        there is an arrow from v2 to v1.
        All these variables have objective value 0.
        The branching priority of these variables can also be set (default = 0)
        '''
        arrow = {}
        n = 0
        for (v1, v2) in permutations(self._bn_variables, 2):
            v = self.addVar(vtype=GRB.BINARY)
            v.BranchPriority = branch_priority
            arrow[v1, v2] = v
            n += 1
        if args.verbose:
            print('%d arrow variables declared' % n, file=sys.stderr)
        self._arrow = arrow

    def add_variables_total_order(self, branch_priority=0):
        '''Adds binary "total order" variables to the model
        and creates a dictionary mapping each ordered pair of BN variables to its
        corresponding total order variable.
        The total order variable corresponding to (v1,v2) is set to 1 iff
        v2 < v1 in the total order.
        All these variables have objective value 0.
        The branching priority of these variables can also be set (default = 0)
        '''
        total_order = {}
        n = 0
        for (v1, v2) in permutations(self._bn_variables, 2):
            v = self.addVar(vtype=GRB.BINARY)
            v.BranchPriority = branch_priority
            total_order[v1, v2] = v
            n += 1
        if args.verbose:
            print('%d total_order variables declared' % n, file=sys.stderr)
        self._total_order = total_order

    def add_variables_adj(self, branch_priority=0):
        '''Adds binary (undirected) "adj" variables to the model
        and creates a dictionary mapping each pair of BN variables to its
        corresponding adj variable.
        The adj variable corresponding to {v1,v2} is set to 1 iff
        there is an arrow from v1 to v2 or an arrow from v2 to v1.
        All these variables have objective value 0.
        The branching priority of these variables can also be set (default = 0)
        '''
        adj = {}
        n = 0
        for (v1, v2) in combinations(self._bn_variables, 2):
            v = self.addVar(vtype=GRB.BINARY)
            v.BranchPriority = branch_priority
            adj[frozenset([v1, v2])] = v
            n += 1
        if args.verbose:
            print('%d adj variables declared' % n, file=sys.stderr)
        self._adj = adj

    def add_variables_gen(self, branch_priority=0):
        '''Adds integer "generation" variables to the model
        and creates a dictionary mapping each BN variable to its
        corresponding generation variable.
        All these variables have objective value 0.
        The branching priority of these variables can also be set (default = 0)
        '''
        gen = {}
        n = 0
        max_num_before = self._n - 1
        for bn_variable in self._bn_variables:
            v = self.addVar(
                ub=max_num_before,
                lb=0,
                vtype=GRB.INTEGER
            )
            gen[bn_variable] = v
            v.BranchPriority = branch_priority
            n += 1
        if args.verbose:
            print('%d generation variables declared' % n, file=sys.stderr)
        self._gen = gen

    def add_variables_gendiff(self, branch_priority=0):
        '''Adds variables representing the difference in generation number
        between distinct BN variables
        '''
        gendiff = {}
        max_num_before = self._n - 1
        n = 0
        for i, v1 in enumerate(self._bn_variables):
            for v2 in self._bn_variables[i + 1:]:
                v = self.addVar(
                    ub=max_num_before,
                    lb=-max_num_before,
                    vtype=GRB.INTEGER)
                v.BranchPriority = branch_priority
                gendiff[v1, v2] = v
                n += 1
        if args.verbose:
            print('%d gen difference variables declared' % n, file=sys.stderr)
        self._gendiff = gendiff

    def add_variables_absgendiff(self, branch_priority=0):
        '''Adds variables representing the absolute difference in generation number
        between distinct BN variables
        '''
        absgendiff = {}
        max_num_before = self._n - 1
        n = 0
        for (v1, v2) in combinations(self._bn_variables, 2):
            v = self.addVar(
                ub=max_num_before,
                lb=1,
                vtype=GRB.INTEGER)
            v.BranchPriority = branch_priority
            absgendiff[frozenset([v1, v2])] = v
            n += 1
        if args.verbose:
            print('%d absolute gen difference variables declared' % n, file=sys.stderr)
        self._absgendiff = absgendiff

    def add_basic_variables(self):
        '''Adds the most useful IP variables
        '''
        self.add_variables_family()
        # self.add_variables_family(best_first_branch_priority=True)
        self.add_variables_arrow(branch_priority=10)
        self.add_variables_adj(branch_priority=10)
        # self.add_variables_total_order()
        # self.add_variables_gen()
        # self.add_variables_gendiff()
        # self.add_variables_absgendiff()

    def add_constraints_oneparentset(self):
        '''Adds the constraint that each child has exactly one parent set
        '''
        n = 0
        for parentset_dict in self._family.values():
            self.addConstr(LinExpr([1.0] * len(parentset_dict), parentset_dict.values()),
                           GRB.EQUAL,
                           1)
            n += 1
        if args.verbose:
            print('%d constraints insisting on exactly one parent set for each variable' % n, file=sys.stderr)

    def _extend(self, ss):
        '''private method used by add_constraints_setpacking
        '''
        new_ss = []
        for si, s in ss:
            # s is a tuple of variable indices
            for i in range(si[-1] + 1, self._n):
                ok1 = True
                new_elt = self._bn_variables[i]
                new_s = s | frozenset([new_elt])
                for child in new_s:
                    others = new_s - frozenset([child])
                    ok2 = False
                    for parentset in self._family[child]:
                        if others <= parentset:
                            ok2 = True
                            break
                    if not ok2:
                        ok1 = False
                        break
                if ok1:
                    new_ss.append((si + (i,), new_s))
        return tuple(new_ss)

    def _spc(self, ss):
        '''private method used by add_constraints_setpacking
        '''
        for x, s in ss:
            fvs = []
            # nonfvs = []
            for child in s:
                others = s - frozenset([child])
                for parentset, fv in self._family[child].items():
                    if others <= parentset:
                        fvs.append(fv)
                    # else:
                    #    nonfvs.append(fv)
            # v = self.addVar(vtype=GRB.BINARY)
            # if len(s) == 3:
            #    v.BranchPriority = 10
            # if len(s) == 4:
            #    v.BranchPriority = 10
            # self.addConstr(LinExpr([1]*len(fvs),fvs),GRB.EQUAL,v)
            self.addConstr(LinExpr([1] * len(fvs), fvs), GRB.LESS_EQUAL, 1)
            # self.addConstr(LinExpr([1]*len(nonfvs),nonfvs),GRB.GREATER_EQUAL,len(s)-1)

    def add_constraints_setpacking(self):
        '''Adds constraints like a<-b,c + b<-a,c + c<-a,b <= 1
        That is an example for a "triple". Also add similar constraints for 4-tuples
        '''
        singles = []
        for i, v in enumerate(self._bn_variables):
            singles.append(((i,), frozenset([v])))
        pairs = self._extend(singles)
        triples = self._extend(pairs)
        quads = self._extend(triples)
        self._spc(triples)
        self._spc(quads)
        if args.verbose:
            print('%d set packing constraints declared' % (len(triples) + len(quads)), file=sys.stderr)

    def add_constraints_arrow_family(self):
        '''Adds constraints linking arrows to family variables
        '''
        n = 0
        m = 0
        for (v1, v2), arrow_var in list(self._arrow.items()):
            vs = [arrow_var]
            non_vs = [arrow_var]
            vals = [-1]
            non_vals = [1]
            for parentset, fv in self._family[v1].items():
                if v2 in parentset:
                    vs.append(fv)
                    vals.append(1)
                else:
                    non_vs.append(fv)
                    non_vals.append(1)
            if len(vs) == 1:
                # arrow can never occur
                self.remove(arrow_var)
                del self._arrow[v1, v2]
                m += 1
            else:
                self.addConstr(LinExpr(vals, vs), GRB.EQUAL, 0)
                self.addConstr(LinExpr(non_vals, non_vs), GRB.EQUAL, 1)
                n += 2
        if args.verbose:
            print('%d constraints linking arrows to family variables declared' % n, file=sys.stderr)
            print('%d arrow variables removed' % m, file=sys.stderr)

    def add_constraints_arrow_total_order(self):
        '''Adds constraints linking arrows to total order variables
        '''
        n = 0
        for (v1, v2), arrow_var in self._arrow.items():
            self.addConstr(arrow_var <= self._total_order[v1, v2])
            n += 1
        if args.verbose:
            print('%d constraints linking arrows to total order variables declared' % n, file=sys.stderr)

    def add_constraints_total_order(self, lazy=False):
        '''Adds constraints so that total order variables
        represent a total order
        '''
        n = 0
        for (v1, v2) in combinations(self._bn_variables, 2):
            self.addConstr(self._total_order[v1, v2] + self._total_order[v2, v1] == 1)
            n += 1
        if lazy:
            self._total_order_lazy = True
        else:
            self._total_order_lazy = False
            for v1, v2 in permutations(self._bn_variables, 2):
                for v3 in self._bn_variables:
                    if v3 != v1 and v3 != v2:
                        self.addConstr(
                            self._total_order[v1, v2] + self._total_order[v2, v3] + self._total_order[v3, v1] <= 2)
                        n += 1
        if args.verbose:
            print('%d total order constraints declared' % n, file=sys.stderr)

    def add_constraints_arrow_adj(self):
        '''Add constraints that there is an undirected adj between
        v1 and v2 in the undirected skeleton if there is either an 
        arrow from v1 to v2, or an arrow from v2 to v1 in the DAG
        '''
        n = 0
        m = 0
        for pair, adj_var in list(self._adj.items()):
            [v1, v2] = list(pair)
            try:
                self.addConstr(
                    LinExpr([-1, 1, 1],
                            [adj_var, self._arrow[v1, v2], self._arrow[v2, v1]]),
                    GRB.EQUAL, 0)
                n += 1
            except KeyError:
                self.remove(adj_var)
                del self._adj[pair]
                m += 1
        if args.verbose:
            print('%d constraints linking arrows to adjs declared' % n, file=sys.stderr)
            print('%d adjacency variables removed' % m, file=sys.stderr)

    def add_constraints_sumgen(self):
        '''Adds the constraint that sum of gen numbers is n*(n-1)/2
        '''
        n = self._n
        self.addConstr(
            LinExpr([1.0] * n, self._gen.values()),
            GRB.EQUAL,
            n * (n - 1) / 2
        )
        if args.verbose:
            print('1 constraint that the sum of gen numbers is', n * (n - 1) / 2, file=sys.stderr)

    def add_constraints_gen_arrow_indicator(self):
        '''Adds constraint stating that an arrow from parent to child
        means that the child's gen number is strictly greater than the parent's
        '''
        n = 0
        for (v1, v2), arrow_var in self._arrow.items():
            self.addConstr((arrow_var == 1) >> (self._gen[v1] <= self._gen[v2] - 1))
            n += 1
        if args.verbose:
            print('%d indicator constraints linking arrows to gen numbers  declared' % n, file=sys.stderr)

    def add_constraints_gendiff(self):
        '''Adds constraints linking gen diff variables to gen variables
        '''
        n = 0
        for (v1, v2), gendiffvar in self._gendiff.items():
            self.addConstr(
                LinExpr([1, -1, -1], [self._gen[v1], self._gen[v2], gendiffvar]),
                GRB.EQUAL,
                0)
            n += 1
        if args.verbose:
            print('%d constraints linking gen diff to gen declared' % n, file=sys.stderr)

    def add_constraints_absgendiff(self):
        '''Adds constraints linking gen diff variables to absolute gen diff variables
        '''
        n = 0
        for (v1, v2), gendiffvar in self._gendiff.items():
            self.addConstr(self._absgendiff[frozenset([v1, v2])] == abs_(gendiffvar))
            n += 1
        if args.verbose:
            print('%d indicator constraints linking abs gen diff to gen diff variables  declared' % n, file=sys.stderr)

    def add_constraints_chordal(self):
        '''Simple constraints to rule out non-chordal DAGs
        i.e. those without v-structures (aka immoralities)
        '''
        n = 0
        dkt = {}
        for child, parent_dict in self._family.items():
            for parent_set, fv in parent_dict.items():
                for pa1, pa2 in combinations(parent_set, 2):
                    ci = frozenset([child, pa1, pa2])
                    try:
                        dkt[ci].append(fv)
                    except KeyError:
                        dkt[ci] = [fv]
        for ci, fvlist in dkt.items():
            for pair in combinations(ci, 2):
                self.addConstr(LinExpr([1] * len(fvlist), fvlist), GRB.LESS_EQUAL, self._adj[frozenset(pair)])
                n += 1
        if args.verbose:
            print('%d constraints ruling out immoralities declared' % n, file=sys.stderr)

    def add_constraints_clusters(self, cluster_cuts=True, matroid_cuts=False, matroid_constraints=False):
        '''Allows cluster constraints to be (lazily) added. If cluster_cuts is
        True they are also added as cuts. If matroid_cuts is true cuts corresponding
        to rank 2 matroid also added
        '''
        self._enforcing_cluster_constraints = True
        self._enforcing_matroid_constraints = matroid_constraints
        self._adding_cluster_cuts = cluster_cuts
        self._adding_matroid_cuts = matroid_cuts
        self.Params.LazyConstraints = 1
        print('(Lazy) "cluster" constraints in use', file=sys.stderr)

    def add_basic_constraints(self):
        '''Add the most useful constraints
        '''
        self.add_constraints_oneparentset()
        self.add_constraints_setpacking()
        self.add_constraints_arrow_family()
        self.add_constraints_arrow_adj()
        # self.add_constraints_arrow_total_order()
        # self.add_constraints_total_order()
        self.add_constraints_clusters()
        # self.add_constraints_sumgen()
        # self.add_constraints_gen_arrow_indicator()
        # self.add_constraints_gendiff()
        # self.add_constraints_absgendiff()

    def mycallback(self, where):
        '''callback for adding cuts and lazy constraints
        '''
        max_cluster_size = self._max_cluster_size
        if (where == GRB.Callback.MIPNODE and
                self.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL and
                self.cbGet(GRB.Callback.MIPNODE_NODCNT) == 0):

            # optionally don't look for cuts if this has already been done
            # 'too' often
            if (self._user_cuts_rounds_limit is not None and
                    self._user_cuts_rounds_count > self._user_cuts_rounds_limit):
                return

            # optionally don't look for cuts if the current bound is not sufficiently
            # distant from the bound we had last time we looked for cuts
            # this is to prevent stalling
            if self._user_cuts_stalling_limit is not None:
                current_bound = self.cbGet(GRB.Callback.MIPNODE_OBJBND)
                if self._last_bound - current_bound < self._user_cuts_stalling_limit:
                    return
                self._last_bound = current_bound

            if self._adding_cluster_cuts:
                self.subip(cutting=True, max_cluster_size=max_cluster_size)
                self._user_cuts_rounds_count += 1
            if self._adding_matroid_cuts:
                self.matroid_subip(cutting=True)
        elif where == GRB.Callback.MIPSOL:

            # optionally don't look for constraints if this has already been done
            # 'too' often - used for solving relaxed versions of the problem
            if (self._user_enforcement_rounds_limit is not None and
                    self._user_enforcement_rounds_count > self._user_enforcement_rounds_limit):
                return

            is_a_dag = True
            if self._enforcing_cluster_constraints:
                is_a_dag = not self.subip(cutting=False, max_cluster_size=max_cluster_size)
                self._user_enforcement_rounds_count += 1
            if self._total_order_lazy == True:
                for v1, v2 in permutations(self._bn_variables, 2):
                    if self.cbGetSolution(self._total_order[v1, v2]) > 0.5:
                        for v3 in self._bn_variables:
                            if v3 == v1 or v3 == v2:
                                continue
                            if (self.cbGetSolution(self._total_order[v2, v3]) > 0.5
                                    and self.cbGetSolution(self._total_order[v3, v1]) > 0.5):
                                self.cbLazy(self._total_order[v1, v2] + self._total_order[v2, v3] + self._total_order[
                                    v3, v1] <= 2)
            if not is_a_dag and self._enforcing_matroid_constraints:
                self.matroid_subip(cutting=False)
            if is_a_dag and self._one_dag_per_MEC == True:
                self.enforce_MEC_representative()

    def matroid_subip(self, cutting):
        '''experimental: generate rank 2 matroid cuts but too slow
        and cuts do not seem too useful
        '''
        matroid_subip = Model("matroid subip")
        matroid_subip.Params.OutputFlag = 0
        matroid_subip.ModelSense = -1
        matroid_subip.Params.Cutoff = -1.5
        matroid_subip.Params.PoolSolutions = 200
        # matroid_subip.Params.PoolSearchMode = 1
        matroid_subip.Params.TimeLimit = 5
        matroid_subip.Params.MIPFocus = 1

        # Dictionaries for matroid subIP variables
        y = {}  # ground set indicators
        sub_fvs = {}  # indicates whether a family variable is in the cut
        circuits2 = {}  # indicates whether a subset of size 2 is a circuit
        not_circuits2 = {}  # indicates whether a subset of size 2 is not a circuit

        try:
            bn_variables_pairs = self._bn_variables_pairs
            bn_variables_pairs_set = self._bn_variables_pairs_set
            bn_variables_triples_of_pairs = self._bn_variables_triples_of_pairs
        except AttributeError:
            bn_variables_pairs = tuple([x for x in combinations(self._bn_variables, 2)])
            bn_variables_pairs_set = tuple([frozenset(x) for x in bn_variables_pairs])
            bn_variables_triples_of_pairs = []
            for (v1, v2, v3) in combinations(self._bn_variables, 3):
                bn_variables_triples_of_pairs.append((frozenset([v1, v2]), frozenset([v1, v3]), frozenset([v2, v3])))
            bn_variables_triples_of_pairs = tuple(bn_variables_triples_of_pairs)
            self._bn_variables_pairs = bn_variables_pairs
            self._bn_variables_pairs_set = bn_variables_pairs_set
            self._bn_variables_triples_of_pairs = bn_variables_triples_of_pairs

        for v in self._bn_variables:
            y[v] = matroid_subip.addVar(vtype=GRB.BINARY, obj=-1)

        for pairset in self._bn_variables_pairs_set:
            circuits2[pairset] = matroid_subip.addVar(vtype=GRB.BINARY)
            not_circuits2[pairset] = matroid_subip.addVar(vtype=GRB.BINARY)

        if cutting:
            fv_vals = self.cbGetNodeRel(self._fv_list)
        else:
            fv_vals = self.cbGetSolution(self._fv_list)

        for i, relval in enumerate(fv_vals):
            if relval > 0:
                child = self._child_list[i]
                parentset = self._parentset_list[i]
                # print child, '<-', ','.join(sorted(parentset)), relval
                v = matroid_subip.addVar(vtype=GRB.BINARY, obj=relval)
                try:
                    sub_fvs[child][parentset] = v
                except KeyError:
                    sub_fvs[child] = {parentset: v}

        matroid_subip.update()

        matroid_subip.addConstr(LinExpr([1] * len(y), y.values()), GRB.GREATER_EQUAL, 4)
        # matroid_subip.addConstr(LinExpr([1]*len(y),y.values()), GRB.LESS_EQUAL, 6)

        for i, (v1, v2) in enumerate(self._bn_variables_pairs):
            pairset = self._bn_variables_pairs_set[i]
            matroid_subip.addConstr(circuits2[pairset] + not_circuits2[pairset] <= y[v1])
            matroid_subip.addConstr(circuits2[pairset] + not_circuits2[pairset] <= y[v2])
            matroid_subip.addConstr(y[v1] + y[v2] <= 1 + circuits2[pairset] + not_circuits2[pairset])

        for pairset1, pairset2, pairset3 in self._bn_variables_triples_of_pairs:
            matroid_subip.addConstr(circuits2[pairset1] + circuits2[pairset2] + not_circuits2[pairset3] <= 2)
            matroid_subip.addConstr(circuits2[pairset1] + circuits2[pairset3] + not_circuits2[pairset2] <= 2)
            matroid_subip.addConstr(circuits2[pairset2] + circuits2[pairset3] + not_circuits2[pairset1] <= 2)

        # want at least one circuit, call it {a,b,c} of size 3, so none of {a,b}, {a,c} and {b,c} can be circuits
        # the following constraint is not sufficient to ensure the existence of a circuit of size 3
        matroid_subip.addConstr(3, GRB.LESS_EQUAL, LinExpr([1] * len(not_circuits2), not_circuits2.values()))

        for child, parent_dict in sub_fvs.items():
            for parent_set, fv in parent_dict.items():
                matroid_subip.addConstr(fv, GRB.LESS_EQUAL, y[child])
                n = len(parent_set)
                matroid_subip.addConstr(
                    fv, GRB.LESS_EQUAL,
                    LinExpr([1] * (n + n * (n - 1) / 2),
                            [circuits2[frozenset([child, parent])] for parent in parent_set] +
                            [not_circuits2[frozenset([pa1, pa2])] for (pa1, pa2) in combinations(parent_set, 2)])
                )

        matroid_subip.optimize()

        # if no matroid constraint found ...
        if matroid_subip.Status == GRB.CUTOFF:
            # print 'no matroid cuts found'
            return False

        for i in range(matroid_subip.Solcount):
            matroid_subip.Params.SolutionNumber = i

            ground_list = [yi for (yi, yv) in y.items() if yv.Xn > 0.5]
            size2_circuits = [circuit for (circuit, v) in circuits2.items() if v.Xn > 0.5]
            size3_circuits = []
            for triple in combinations(ground_list, 3):
                triple = frozenset(triple)
                for c in size2_circuits:
                    if c <= triple:
                        break
                else:
                    size3_circuits.append(triple)

            # if ended up without any circuits of size 3
            # then no good
            if size3_circuits == []:
                continue

            size23_circuits = size2_circuits + size3_circuits

            # check that the matroid is connected
            connected = True
            for (v1, v2) in combinations(ground_list, 2):
                pair = frozenset((v1, v2))
                for circuit in size23_circuits:
                    if pair <= circuit:
                        break
                else:
                    connected = False
                    break

            if not connected:
                continue

            # print
            # print 'Ground set', ground_list
            # print 'Circuits 2', size2_circuits
            # print 'Circuits 3', size3_circuits
            # print matroid_subip.objVal

            rhs = len(ground_list) - 2
            lexpr = LinExpr()
            lhs_string = ''
            activity = 0.0
            for child in ground_list:
                intersecting_fvs = []
                non_intersecting_fvs = []
                for parent_set, fv in self._family[child].items():
                    tmp = parent_set | frozenset([child])
                    for circuit in size23_circuits:
                        if child in circuit and circuit <= tmp:
                            intersecting_fvs.append(fv)
                            lhs_string += ' + {0}<-{1}'.format(child, ','.join(sorted(parent_set)))
                            try:
                                activity += sub_fvs[child][parent_set].obj
                            except KeyError:
                                pass
                            break
                    else:
                        non_intersecting_fvs.append(fv)
                len_i = len(intersecting_fvs)
                len_ni = len(non_intersecting_fvs)
                if len_i > len_ni:
                    lexpr.addTerms([-1] * len_ni, non_intersecting_fvs)
                    rhs -= 1
                else:
                    lexpr.addTerms([1] * len_i, intersecting_fvs)
            if cutting:
                self.cbCut(lexpr, GRB.LESS_EQUAL, rhs)
            else:
                self.cbLazy(lexpr, GRB.LESS_EQUAL, rhs)
            # print lhs_string, '<=', len(ground_list)-2, 'activity = ', activity
        return True

    def vanilla_cons(self, cutting, cluster):
        '''in each cluster one variable must come last
        and thus have all other cluster members available as
        parents
        '''
        lexpr = LinExpr()
        rhs = 1
        for child in cluster:
            child_dkt = self._family[child]
            n = len(self._ordered_parentsets[child])
            for i, parent_set in enumerate(self._ordered_parentsets[child]):
                if parent_set <= cluster:
                    if i < n / 2:
                        lexpr.addTerms([1] * (i + 1), [child_dkt[ps] for ps in self._ordered_parentsets[child][:i + 1]])
                    else:
                        lexpr.addTerms([-1] * (n - (i + 1)),
                                       [child_dkt[ps] for ps in self._ordered_parentsets[child][i + 1:]])
                        rhs -= 1
                    # don't consider parent sets worse than this one
                    break
        if cutting:
            self.cbCut(lexpr, GRB.GREATER_EQUAL, rhs)
        else:
            self.cbLazy(lexpr, GRB.GREATER_EQUAL, rhs)

    def subip(self, cutting, max_cluster_size=None):
        '''Sub-IP for finding cluster cuts which separate the solution
        to the current linear relaxation.
        Returns true iff an efficacious cut/constraint was found
        '''
        # get vals of all family variables in LP relaxation solution
        if cutting:
            fv_vals = self.cbGetNodeRel(self._fv_list)
        else:
            fv_vals = self.cbGetSolution(self._fv_list)
        subip = Model("subip")
        subip.Params.OutputFlag = 0
        subip.ModelSense = -1
        subip.Params.PoolSolutions = 200
        subip.Params.PoolSearchMode = 1
        if cutting:
            subip.Params.TimeLimit = self._subip_cutting_timelimit
        # need to set this strictly above -1 to work
        # in contradiction to gurobi documentation
        subip.Params.Cutoff = self._subip_cutoff
        y = {}
        sub_fvs = {}
        for v in self._bn_variables:
            y[v] = subip.addVar(vtype=GRB.BINARY, obj=-1)
            sub_fvs[v] = {}
        for i, relval in enumerate(fv_vals):
            if relval > 0:
                child = self._child_list[i]
                parentset = self._parentset_list[i]
                # print child, parentset
                v = subip.addVar(vtype=GRB.BINARY, obj=relval)
                try:
                    sub_fvs[child][parentset] = v
                except KeyError:
                    sub_fvs[child] = {parentset: v}
        subip.update()
        if max_cluster_size is not None:
            subip.addConstr(LinExpr([1] * len(y), y.values()), GRB.LESS_EQUAL, max_cluster_size)
        subip.addConstr(LinExpr([1] * len(y), y.values()), GRB.GREATER_EQUAL, 2)
        for child, parent_dict in sub_fvs.items():
            for parentset, fv in parent_dict.items():
                subip.addConstr(fv, GRB.LESS_EQUAL, y[child])
                subip.addConstr(fv, GRB.LESS_EQUAL, LinExpr([(1, y[parent]) for parent in parentset]))
        subip.optimize()

        # if no cluster constraint found ...
        if subip.Status == GRB.CUTOFF:
            return False

        if cutting:
            # add all found cuts
            nsols = subip.Solcount
        else:
            # only add one constraint
            nsols = 1
        for i in range(nsols):
            subip.Params.SolutionNumber = i
            cluster = []
            for v, yv in y.items():
                if yv.Xn > 0.5:
                    cluster.append(v)
            # print 'Cluster', cluster
            cluster_set = frozenset(cluster)
            rhs = len(cluster) - 1
            lexpr = LinExpr()
            for child in cluster:
                intersecting_fvs = []
                non_intersecting_fvs = []
                for parent_set, fv in self._family[child].items():
                    if parent_set & cluster_set:
                        intersecting_fvs.append(fv)
                    else:
                        non_intersecting_fvs.append(fv)
                len_i = len(intersecting_fvs)
                len_ni = len(non_intersecting_fvs)
                # there is exactly one parent set per child
                # use this to minimise number of variables in cut/constraint
                if len_i > len_ni:
                    lexpr.addTerms([-1] * len_ni, non_intersecting_fvs)
                    rhs -= 1
                else:
                    lexpr.addTerms([1] * len_i, intersecting_fvs)
            if cutting:
                self.cbCut(lexpr, GRB.LESS_EQUAL, rhs)
            else:
                self.cbLazy(lexpr, GRB.LESS_EQUAL, rhs)
            # self.vanilla_cons(cutting,frozenset(cluster))
        return True


if __name__ == '__main__':
    model = Gobnilp("gobnilp")
    argv1 = args.scores
    if argv1 == '-':
        # if reading from C gobnilp output, currently have to filter out first 10 lines like this:
        # bin/gobnilp -x data/alarm_1000.dat | tail -n +10 | gurobi.sh simple_gobnilp.py -
        fobj = sys.stdin
    else:
        fobj = open(argv1)
    model.read_local_scores(fobj)
    model.add_basic_variables()
    #model.add_variables_local_scores()
    model.update()
    model.Params.PoolSolutions = args.nsols  # save k solutions
    #model.Params.BranchDir = 1
    if args.kbest:
        model.Params.PoolSearchMode = 2  # find k best solutions
    model.add_basic_constraints()
    if args.mec:
        model.add_constraints_one_dag_per_MEC()
    # if args.chordal:
    #    model.add_constraints_chordal()
    #    model._max_cluster_size = 3
    model.optimize(lambda model, where: model.mycallback(where))
    for i in range(model.Solcount):
        model.Params.SolutionNumber = i
        model.print_simple_output()
