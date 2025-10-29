import polars as PL
from ortools.sat.python import cp_model
from ortools.graph.python import linear_sum_assignment
"""
version of the package
"""
version='1.0'
"""
date of this version of the package
"""
date='2025-10-29'
class multipletsError(Exception):
  """Exception raised by package multiplets"""
  def __init__(self,method:str,text:str):
    """
    Initialize the exception.
    
    Parameters
    ----------
    method
        Name of the method which raised the exception.
    text
        Text of the exception.
    """
    self.method=method
    self.text=text
  def __str__(self)->str:
    return f'{self.method}: {self.text}'
class multipletsErrorNoSolution(multipletsError):
  """Exception raised by package multiplets when no solution is found"""
  def __init__(self,text:str):
    """
    Initialize the exception.
    
    Parameter
    ---------
    text
        Text of the exception.
    """
    super().__init__('find_multiplets',text)
class _multipletsWarning():
  """Warning issued by package multiplets"""
  def __init__(self,method:str,text:str):
    """
    Initialize and print the warning.
    
    Parameters
    ----------
    method
        Name of the method which issued the warning.
    text
        Text of the warning.
    """
    self.method=method
    self.text=text
    print(f'{self.__class__.__name__}: {self}')
  def __str__(self)->str:
    return f'{self.method}: {self.text}'
class multiplets:
  """The main class for generating a set of multiplets"""
  def __init__(self,df:PL.DataFrame,id:str,group:str,*,check_too_big:bool=True,colname_rowid:str='_rowid'):
    """
    Initialize the data structures. Partition vertices into groups.
    
    Parameters
    ----------
    df
        Dataframe containing columns id and group where each row represents one person. Dataframe can contain another columns and they can be used when generating the multiplets.
    id
        Name of the column containing the personal identifiers. Values in this column must be nonempty and unique.
    group
        Name of the column which contains values, according to which the set of persons is split into groups.
    check_too_big
        If check_too_big is True and the column group has more than 10 values then the constructor raises an error. This can happen when there is some user mistake.
    colname_rowid
        Name of a column which is used for internal manipulation with dataframes.
    """
    self.init_OK=False
    if df.height==0:
      raise multipletsError('__init__','Dataframe has zero rows! Exiting.')
    if id in df.columns:
      self.colname_id=id
    else:
      raise multipletsError('__init__',f'Dataframe does not contain column {id}! Exiting.')
    if group in df.columns:
      self.colname_group=group
    else:
      raise multipletsError('__init__',f'Dataframe does not contain column {group}! Exiting.')
    if df.select(id).height!=df.select(id).unique().drop_nulls().height:
      raise multipletsError('__init__',f'Values in column {id} are not nonemply and unique! Exiting.')
    self.colname_rowid=colname_rowid
    self.vertices=df.partition_by(group,include_key=False,as_dict=True)
    t=list(self.vertices.keys())
    if (None,) in t:
      t.remove((None,))
      self.part_keys=sorted(t)+[(None,)]
    else:
      self.part_keys=sorted(t)
    if len(self.part_keys)==1:
      raise multipletsError('__init__',f'Partition by {group} has only one part! Exiting.')
    elif (len(self.part_keys)>10) and check_too_big:
      raise multipletsError('__init__',f'Partition by {group} has too much ({len(self.part_keys)}) parts! Exiting.')
    if len(self.part_keys)==2:
      for k in self.vertices:
        self.vertices[k]=self.vertices[k].with_row_index(self.colname_rowid)
    self._was_col_warning=False
    self.init_OK=True
  def _col(self,strorexpr:'Expr|str|int|float')->'Expr':
    """
    Convert argument into Polars expression, depending on the type of the argument.
    
    Parameter
    ---------
    strorexpr
        type Expr - return the argument
        type str - return column with this name
        type int or float - return column with this literal value
        otherwise - return column with literal 0, issuing a warning
    """
    if isinstance(strorexpr,PL.expr.expr.Expr):
      return strorexpr
    elif isinstance(strorexpr,str):
      return PL.col(strorexpr)
    elif isinstance(strorexpr,(int,float)):
      return PL.lit(strorexpr)
    else:
      if not self._was_col_warning:
        _multipletsWarning('_col',f'Argument {strorexpr} has unsupported type {type(strorexpr)}! Using 0 instead.')
        self._was_col_warning=True
      return PL.lit(0)
  def init_edges(self,weight:'Expr|str|int|float'=0,filter:'Expr|int|float|NoneType'=None,*,agg_horizontal:'function'=PL.sum_horizontal,colname_weight:str='_weight')->PL.DataFrame:
    """
    Initialize edges and hyperedges between vertices.
    
    An undirected edge connects persons from different groups which satisfy the condition of parameter filter. The weight of an edge is given by parameter weight. The edges are stored as a dataframe in atribute edges.
    
    A hyperedge connects only persons across all groups such that each pair of these persons is connected by an edge. The weight of a hyperedge is computed by function agg_horizontal applied to weights of all these edges.
    
    The method returns a dataframe where each row represents one hyperedge, with column '_weight' representing a weight of this hyperedge. It is ensured that for the same input this dataframe will be the same. The returned dataframe is also stored in atribute hyperedges.
    
    Parameters
    ----------
    weight
        The weight of an edge. It can be name of a column or expression made by other columns. If the original dataframe in multiplets.__init__ contains another columns then the names of these columns are suffixed by '_A' and '_B', respectively  for the two persons connected by the edge. These suffixed column names can also be used in the expression.
        The value of weight is temporarily stored in column '_weight'.
    filter
        A condition which must be satisfied by all edges. If filter is some number x of type int or float then the condition is
            pl.col('_weight')<=x
        Implicitely None which means that no filter is applied.
    agg_horizontal
        An aggregation function used to compute weight of a hyperedge. This function is applied to weights of all edges contained in the hyperedge.
    colname_weight
        Name of a column which is used for internal manipulation with dataframes.
    
    Example
    -------
    Suppose that the persons are in dataframe
           ┌──┬─────┬─────┐
           │id┆group┆value│
        df=╞══╪═════╪═════╡
           │ 1┆ 'D' ┆  10 │
           │ 2┆ 'D' ┆  20 │
           │ 3┆ 'E' ┆  30 │
           │ 4┆ 'E' ┆  40 │
           │ 5┆ 'F' ┆  50 │
           │ 6┆ 'F' ┆  60 │
           └──┴─────┴─────┘
    Then (assuming that package polars was imported as pl)
        mp=multiplets(df,'id','group')
        mp.init_edges(weight=(pl.col('value_A')-pl.col('value_B')).abs(), filter=30)
    returns
                      ┌────┬────┬────┬───────┐
                      │id_0┆id_1┆id_2┆_weight│
        mp.hyperedges=╞════╪════╪════╪═══════╡
                      │  2 ┆  3 ┆  5 ┆   60  │
                      │  2 ┆  4 ┆  5 ┆   60  │
                      └────┴────┴────┴───────┘
    Here in the first row we have abs(20-30)+abs(20-50)+abs(30-50)=60.
    On the other hand,
        mp=multiplets(df,'id','group')
        mp.init_edges(weight=(pl.col('value_A')-pl.col('value_B')).abs(), filter=30,
            agg_horizontal=pl.max_horizontal)
    returns
                      ┌────┬────┬────┬───────┐
                      │id_0┆id_1┆id_2┆_weight│
        mp.hyperedges=╞════╪════╪════╪═══════╡
                      │  2 ┆  3 ┆  5 ┆   30  │
                      │  2 ┆  4 ┆  5 ┆   30  │
                      └────┴────┴────┴───────┘
    Here in the first row we have max(abs(20-30),abs(20-50),abs(30-50))=30.
    
    Example
    -------
    In the case that there are only two groups, dataframe mp.hyperedges contains also columns '_rowid_0' and '_rowid_1' which are used later for faster search.
    Suppose that the persons are in dataframe
           ┌──┬─────┬─────┐
           │id┆group┆value│
        df=╞══╪═════╪═════╡
           │ 1┆ 'D' ┆  10 │
           │ 2┆ 'D' ┆  20 │
           │ 3┆ 'E' ┆  30 │
           │ 4┆ 'E' ┆  40 │
           └──┴─────┴─────┘
    Then (assuming that package polars was imported as pl)
        mp=multiplets(df,'id','group')
        mp.init_edges(weight=(pl.col('value_A')-pl.col('value_B')).abs(), filter=20)
    returns
                      ┌────────┬────────┬────┬────┬───────┐
                      │_rowid_0┆_rowid_1┆id_0┆id_1┆_weight│
        mp.hyperedges=╞════════╪════════╪════╪════╪═══════╡
                      │    0   ┆    0   ┆  1 ┆  3 ┆   20  │
                      │    1   ┆    0   ┆  2 ┆  3 ┆   10  │
                      │    1   ┆    1   ┆  2 ┆  4 ┆   20  │
                      └────────┴────────┴────┴────┴───────┘
    Here in the first row we have abs(10-30)=20.
    """
    if not self.init_OK:
      raise multipletsError('init_edges','Instance is not initialized correctly! Exiting.')
    self.colname_weight=colname_weight
    if filter is None:
      filter=True
    elif isinstance(filter,(int,float)):
      filter=(PL.col(self.colname_weight)<=filter)
    self.edges={}
    for i in range(len(self.part_keys)-1):
      for j in range(i+1,len(self.part_keys)):
        self.edges[i,j]=(
          self.vertices[self.part_keys[i]].select(PL.all().name.suffix('_A'))
            .join(self.vertices[self.part_keys[j]].select(PL.all().name.suffix('_B')), how='cross')
            .with_columns(self._col(weight).alias(self.colname_weight))
            .filter(filter)
            .select(
              ([PL.col(f'{self.colname_rowid}_A').alias(f'{self.colname_rowid}_{i}'),
                  PL.col(f'{self.colname_rowid}_B').alias(f'{self.colname_rowid}_{j}')]
                if len(self.part_keys)==2 else [])+
              [PL.col(f'{self.colname_id}_A').alias(f'{self.colname_id}_{i}'),
                PL.col(f'{self.colname_id}_B').alias(f'{self.colname_id}_{j}'),
                PL.col(self.colname_weight).alias(f'{self.colname_weight}_{i}_{j}')]))
    first_join=True
    for i in range(len(self.part_keys)-1):
      for j in range(i+1,len(self.part_keys)):
        if first_join:
          self.hyperedges=self.edges[i,j]
          first_join=False
        else:
          self.hyperedges=self.hyperedges.join(
            self.edges[i,j], on=[f'{self.colname_id}_{i}', f'{self.colname_id}_{j}'][:i+1], how='inner')
    self.hyperedges=(self.hyperedges
      .select(
        ([f'{self.colname_rowid}_{i}' for i in range(len(self.part_keys))] if len(self.part_keys)==2 else [])+
        [f'{self.colname_id}_{i}' for i in range(len(self.part_keys))]+
        [agg_horizontal(f'{self.colname_weight}_{i}_{j}'
            for i in range(len(self.part_keys)-1) for j in range(i+1,len(self.part_keys))
          ).alias(self.colname_weight)])
      .sort([f'{self.colname_id}_{i}' for i in range(len(self.part_keys))]))
    return self.hyperedges
  def join(self,left_df:PL.DataFrame,right_df:PL.DataFrame,right_on:'str|NoneType'=None)->PL.DataFrame:
    """
    Left join a dataframe containing multiplets with horizontal id's, with another dataframe containing id's and other values. In the returned dataframe each row represents a multiplet, the dataframe contains horizontal id's and horizontal values.
    
    Parameters
    ----------
    left_df
        A dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. The id's must be in the form 'id_0', 'id_1' etc. where 'id' is the name of the id column used in the constructor.
    right_df
        A dataframe where each row represents a person, with vertical id's and possibly some other columns.
    right_on
        Name of the column in right_df which contains the id's.
    
    Example
    -------
    Suppose that dataframe df contains in column 'group' three different values. Let
        mp=multiplets(df,'id','group')
    Then
                        ┌────┬────┬────┐           ┌──┬─────┐
                        │id_0┆id_1┆id_2│           │ID┆value│
        mp.join(left_df=╞════╪════╪════╡, right_df=╞══╪═════╡, right_on='ID')
                        │  1 ┆  2 ┆  3 │           │ 1┆  10 │
                        │  4 ┆  5 ┆  6 │           │ . . .  │
                        └────┴────┴────┘           │ 6┆  60 │
                                                   └──┴─────┘
    returns
        ┌────┬────┬────┬───────┬───────┬───────┐
        │id_0┆id_1┆id_2┆value_0┆value_1┆value_2│
        ╞════╪════╪════╪═══════╪═══════╪═══════╡
        │  1 ┆  2 ┆  3 ┆   10  ┆   20  ┆   30  │
        │  4 ┆  5 ┆  6 ┆   40  ┆   50  ┆   60  │
        └────┴────┴────┴───────┴───────┴───────┘
    """
    if not self.init_OK:
      raise multipletsError('join','Instance is not initialized correctly! Exiting.')
    if right_on is None:
      right_on=self.colname_id
    if right_on not in right_df.columns:
      raise multipletsError('join',f'Right dataframe does not contain column {right_on}! Exiting.')
    res=left_df
    for i in range(len(self.part_keys)):
      res=res.join(right_df.rename(lambda s:f'{s}_{i}'),
        left_on=f'{self.colname_id}_{i}', right_on=f'{right_on}_{i}', how='left')
    return res
  def unpivot(self,df:PL.DataFrame,*,all_columns:bool=False)->PL.DataFrame:
    """
    Unpivot a dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. In the returned dataframe each row represents a person, the dataframe contains vertical id's, position in the multiplets (column named 'multiplet_id' and column named in the same way as 'group' in the constructor) and possibly the other columns.
    
    Parameters
    ----------
    df
        A dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. The id's must be in the form 'id_0', 'id_1' etc. where 'id' is the name of the id column used in the constructor.
    all_columns
        If set to True then values in the other columns will be copied fror a multiplet to each its person.
    
    Example
    -------
    Suppose that dataframe df contains in column 'group' unique values 'A', 'B', 'C'. Let
        mp=multiplets(df,'id','group')
    Then
                      ┌────┬────┬────┬─────┐
                      │id_0┆id_1┆id_2┆value│
        mp.unpivot(df=╞════╪════╪════╪═════╡, all_columns=True)
                      │  1 ┆  2 ┆  3 ┆  10 │
                      │  4 ┆  5 ┆  6 ┆  20 │
                      └────┴────┴────┴─────┘
    returns
        ┌────────────┬─────┬──┬─────┐
        │multiplet_id┆group┆id┆value│
        ╞════════════╪═════╪══╪═════╡
        │      0     ┆ 'A' ┆ 1┆  10 │
        │      0     ┆ 'B' ┆ 2┆  10 │
        │      0     ┆ 'C' ┆ 3┆  10 │
        │      1     ┆ 'A' ┆ 4┆  20 │
        │      1     ┆ 'B' ┆ 5┆  20 │
        │      1     ┆ 'C' ┆ 6┆  20 │
        └────────────┴─────┴──┴─────┘
    In case all_columns=False, only columns 'multiplet_id', 'group' and 'id' are kept.
    """
    if not self.init_OK:
      raise multipletsError('unpivot','Instance is not initialized correctly! Exiting.')
    t=df.with_row_index('multiplet_id')
    res=(t.unpivot(
        [f'{self.colname_id}_{i}' for i in range(len(self.part_keys))],
        index='multiplet_id',
        variable_name='_group',value_name=self.colname_id)
      .join(PL.DataFrame({
            '_group': [f'{self.colname_id}_{i}' for i in range(len(self.part_keys))],
            self.colname_group: [self.part_keys[i][0] for i in range(len(self.part_keys))]}),
          on='_group', how='left')
      .select('multiplet_id',self.colname_group,self.colname_id))
    if all_columns:
      res=res.join(
        t.drop([f'{self.colname_id}_{i}' for i in range(len(self.part_keys))]),
        on='multiplet_id', how='left')
    return res
  def find_multiplets(self,df:'PL.DataFrame|NoneType'=None,*,weight='_weight',multiplier:'int|float|NoneType'=None,penalty:'int|float|NoneType'=None,force_CPSAT:bool=False,max_time:'int|float|NoneType'=None,verbose:int=0):
    """
    Try to find the best set of multiplets from a given set of possible multiplets.
    
    The method returns the set with the most hyperedges and in the case of a tie, with the smallest total weight. The returned set is represented by a dataframe where each row represents a multiplet.
    
    Usually there is more than one optimal solution. In this case, each execution of find_multiplets can return a different set of multiplets.
    
    Parameters
    ----------
    df
        A dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. The id's must be in the form 'id_0', 'id_1' etc. where 'id' is the name of the id column used in the constructor.
        The dataframe should contain column given by parameter weight.
        Implicitely df=None which means that self.hyperedges is used.
    weight
        Name of the column, with respect to which the minimization will be done.
        In the case that the dataframe does not contain this column, then this column is filled with zeroes.
    multiplier
        In algorithm LSA, the weights are multiplied by this number before conversion to integers.
        Implicitely multiplier=1.
    penalty
        Penalty added to the total weight in the case that the set of multiplets does not cover the smallest group.
        The penalty is added for every one person from the smallest group which is not covered by the multiplets.
        Implicitely penalty=None which means that the algorithm looks for the maximal number of multiplets, and for this number of multiplets the algorithm looks for the set of multiplets with the smallest total weight.
    force_CPSAT
        For more than two groups, algorithm CPSAT is used for searching.
        Also for two groups and force_CPSAT=True, algorithm CPSAT is used.
        Implicitely force_CPSAT=False which means that for two groups, algorithm LSA is used.
        Algorithm LSA is much faster, but it works only for two groups and for integer values of weights.
    max_time
        Algorithm CPSAT will be interrupted after this number of seconds.
        Implicitely max_time=None which means that the algorithm will not be interrupted.
    verbose
        If verbose>=1 then some detailed information about the solution is printed.
        If verbose>=2 then the information about currently the best sets is printed.
        Implicitely verbose=0.
    """
    if not self.init_OK:
      raise multipletsError('find_multiplets','Instance is not initialized correctly! Exiting.')
    if df is None:
      df=self.hyperedges
    if weight not in df.columns:
      df=df.with_columns(pl.lit(0).alias(weight))
    if penalty is None:
      penalty=int(df.select(weight).max().item()+3)
    if force_CPSAT or (len(self.part_keys)>=3):
      if force_CPSAT and (len(self.part_keys)>=3):
        _multipletsWarning('find_multiplets','Argument force_CPSAT is not used for more than two groups.')
      if multiplier is not None:
        _multipletsWarning('find_multiplets','Argument multiplier is not used for algorithm CPSAT.')
      self.model_type='CPSAT'
      self.model=cp_model.CpModel()
      self.x=[self.model.new_int_var(0,1,f'x{i}') for i in range(df.height)]
      t=self.unpivot(df).group_by(self.colname_id).agg('multiplet_id')
      for i in range(t.height):
        self.model.add(sum(self.x[j] for j in t.item(i,1))<=1)
      self.model.maximize(sum((penalty-w)*self.x[i] for i,w in enumerate(df.get_column(weight))))
      self.solver=cp_model.CpSolver()
      if max_time is not None:
        self.solver.parameters.max_time_in_seconds=max_time
      if verbose>=2:
        self.status=self.solver.solve(self.model,cp_model.ObjectiveSolutionPrinter())
      else:
        self.status=self.solver.solve(self.model)
      if verbose>=1:
        print(f'  status   : {self.solver.status_name(self.status)}')
        print(f'  conflicts: {self.solver.num_conflicts}')
        print(f'  branches : {self.solver.num_branches}')
        print(f'  wall time: {self.solver.wall_time} s')
      if self.status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        if verbose>=1:
          print(f'Maximum of objective function: {self.solver.objective_value}')
        for i in range(len(self.x)):
          if self.solver.value(self.x[i]):
            self.multiplets=PL.concat(df[i] for i in range(len(self.x)) if self.solver.value(self.x[i]))
            break
        else:
          self.multiplets=df.head(0)
          _multipletsWarning('find_multiplets','The solution has zero multiplets.')
      else:
        self.multiplets=df.head(0)
        if verbose>=1:
          print('No solution found!')
        raise multipletsErrorNoSolution('No solution found!')
    else:
      if max_time is not None:
        _multipletsWarning('find_multiplets','Argument max_time is not used for algorithm LSA.')
      self.model_type='LSA'
      if multiplier is None:
        multiplier=1
      heights=[self.vertices[self.part_keys[i]].height for i in range(2)]
      def rangedf(a0,b0,a1,b1,w):
        return (PL.DataFrame({f'{self.colname_rowid}_0':range(a0,b0)})
          .join(PL.DataFrame({f'{self.colname_rowid}_1':range(a1,b1)}), how='cross')
          .with_columns(PL.lit(w).alias(weight)))
      df=(PL.concat(
          [df]+
          ([rangedf(heights[0],heights[1],0,heights[1],0)] if heights[0]<heights[1] else
            [rangedf(0,heights[0],heights[1],heights[0],0)] if heights[1]<heights[0] else [])+
          [rangedf(max(heights),sum(heights),0,heights[1],penalty/2),
            rangedf(0,heights[0],max(heights),sum(heights),penalty/2),
            PL.DataFrame({
              f'{self.colname_rowid}_0':range(max(heights),sum(heights)),
              f'{self.colname_rowid}_1':range(max(heights),sum(heights)),
              weight:0})],
          how='diagonal_relaxed'))
      if verbose>=2:
        PL.Config.set_tbl_rows(100)
        print(f'{df=}')
        PL.Config.set_tbl_rows(10)
      self.assignment=linear_sum_assignment.SimpleLinearSumAssignment()
      self.assignment.add_arcs_with_cost(
        df.get_column(f'{self.colname_rowid}_0'),
        df.get_column(f'{self.colname_rowid}_1'),
        (df.get_column(weight)*multiplier).cast(PL.Int64))
      self.status=self.assignment.solve()
      if self.status==self.assignment.OPTIMAL:
        self.multiplets=(
          PL.DataFrame({f'{self.colname_rowid}_0':range(self.assignment.num_nodes()),
            f'{self.colname_rowid}_1':(self.assignment.right_mate(i) for i in range(self.assignment.num_nodes()))})
          .join(df, on=[f'{self.colname_rowid}_0',f'{self.colname_rowid}_1'], how='left')
          .drop_nulls([f'{self.colname_id}_0',f'{self.colname_id}_1'])
          .select(f'{self.colname_id}_0',f'{self.colname_id}_1',weight))
        if verbose>=2:
          PL.Config.set_tbl_rows(100)
          print(f'{self.multiplets=}')
          PL.Config.set_tbl_rows(10)
        if verbose>=1:
          print(f'Maximum of objective function: {self.assignment.optimal_cost()}')
          print(f'Total weight: {self.multiplets.select(weight).sum().item()}')
      elif self.status==self.assignment.INFEASIBLE:
        self.multiplets=df.head(0)
        if verbose>=1:
          print('No assignment is possible!')
        raise multipletsErrorNoSolution('No assignment is possible!')
      elif self.status==self.assignment.POSSIBLE_OVERFLOW:
        self.multiplets=df.head(0)
        if verbose>=1:
          print('Some input costs are too large and may cause an integer overflow!')
        raise multipletsErrorNoSolution('Some input costs are too large and may cause an integer overflow!')
    return self.multiplets
