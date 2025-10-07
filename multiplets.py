import polars as PL
from ortools.sat.python import cp_model
"""
version of the package
"""
version='0.8'
"""
date of this version of the package
"""
date='2025-10-07'
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
  def __init__(self,df:PL.DataFrame,id:str,group:str,*,check_too_big:bool=True):
    """
    Initialize the data structures. Partition vertices into groups.
    
    Parameters
    ----------
    df
        Dataframe containing columns id and group where each row represents one proband. Dataframe can contain another columns and they can be used when generating the multiplets.
    id
        Name of the column containing the proband identifiers. Values in this column must be unique.
    group
        Name of the colums which contains values, according to which the set of probands is split into groups.
    check_too_big
        If check_too_big is True and the column group has more than 10 values then the constructor raises an error. This can happen when there is some user mistake.
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
    if df.select(id).height!=df.select(id).unique().height:
      raise multipletsError('__init__',f'Values in column {id} are not unique! Exiting.')
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
  def init_edges(self,weight:'Expr|str|int|float'=0,filter:'Expr|int|float'=1,*,agg_horizontal:'function'=PL.sum_horizontal)->PL.DataFrame:
    """
    Initialize edges and hyperedges between vertices.
    
    An undirected edge connects probands from different groups which satisfy the condition of parameter filter. The weight of an edge is given by parameter weight. The edges are stored as a dataframe in atribute edges.
    
    A hyperedge connects only probands across all groups such that each pair of these probands is connected by an edge. The weight of a hyperedge is computed by function agg_horizontal applied to weights of all these edges.
    
    The method returns a dataframe where each row represents one hyperedge, with column '_weight' representing a weight of this hyperedge. It is ensured that for the same input this dataframe will be the same. The returned dataframe is also stored in atribute hyperedges.
    
    Parameters
    ----------
    weight
        The weight of an edge. It can be name of a column or expression made by other columns. If the original dataframe in multiplets.__init__ contains another columns then the names of these columns are suffixed by '_A' and '_B', respectively  for the two probands connected by the edge. These suffixed column names can also be used in the expression.
        The value of weight is temporarily stored in column '_weight'.
    filter
        A condition which must be satisfied by all edges. If filter is some number x of type int or float then the condition is
            pl.col('_weight')<=x
    agg_horizontal
        An aggregation function used to compute weight of a hyperedge. This function is applied to weights of all edges contained in the hyperedge.
    
    Example
    -------
    Suppose that the probands are in dataframe
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
    """
    if not self.init_OK:
      raise multipletsError('init_edges','Instance is not initialized correctly! Exiting.')
    if isinstance(filter,(int,float)):
      filter=(PL.col('_weight')<=filter)
    self.edges={}
    for i in range(len(self.part_keys)-1):
      for j in range(i+1,len(self.part_keys)):
        self.edges[i,j]=(
          self.vertices[self.part_keys[i]].select(PL.all().name.suffix('_A'))
            .join(self.vertices[self.part_keys[j]].select(PL.all().name.suffix('_B')), how='cross')
            .with_columns(self._col(weight).alias('_weight'))
            .filter(filter)
            .select(
              PL.col(f'{self.colname_id}_A').alias(f'{self.colname_id}_{i}'),
              PL.col(f'{self.colname_id}_B').alias(f'{self.colname_id}_{j}'),
              PL.col('_weight').alias(f'_weight_{i}_{j}')))
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
      .select(f'^{self.colname_id}_.*$', agg_horizontal(f'^_weight_.*$').alias('_weight'))
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
        A dataframe where each row represents a proband, with vertical id's and possibly some other columns.
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
    Unpivot a dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. In the returned dataframe each row represents a proband, the dataframe contains vertical id's, position in the multiplets (column named 'multiplet_id' and column named in the same way as 'group' in the constructor) and possibly the other columns.
    
    Parameters
    ----------
    df
        A dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. The id's must be in the form 'id_0', 'id_1' etc. where 'id' is the name of the id column used in the constructor.
    all_columns
        If set to True then values in the other columns will be copied fror a multiplet to each its proband.
    
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
  def find_multiplets(self,df:'PL.DataFrame|NoneType'=None,*,weight='_weight',verbose:int=0):
    """
    Try to find the best set of multiplets from a given set of possible multiplets.
    
    The method returns the set with the most hyperedges and in the case of a tie, with the smallest total weight. The returned set is represented by a dataframe where each row represents a multiplet.
    
    Parameters
    ----------
    df
        A dataframe where each row represents a multiplet, with horizontal id's and possibly some other columns. The id's must be in the form 'id_0', 'id_1' etc. where 'id' is the name of the id column used in the constructor.
        The dataframe should contain column given by parameter weight.
        Implicitely df=None which means that self.hyperedges is used.
    weight
        Name of the column, with respect to which the minimization will be done.
        In the case that the dataframe does not contain this column, then this column is filled with zeroes.
    verbose
        If verbose>0 then the information about currently the best sets is printed.
        Implicitely verbose=0.
    """
    if not self.init_OK:
      raise multipletsError('find_multiplets','Instance is not initialized correctly! Exiting.')
    if df is None:
      df=self.hyperedges
    if weight not in df.columns:
      df=df.with_columns(pl.lit(0).alias(weight))
    self.model=cp_model.CpModel()
    self.x=[self.model.new_int_var(0,1,f'x{i}') for i in range(df.height)]
    t=self.unpivot(df).group_by(self.colname_id).agg('multiplet_id')
    for i in range(t.height):
      self.model.add(sum(self.x[j] for j in t.item(i,1))<=1)
    C=df.select(weight).max().item()*min(d.height for d in self.vertices.values())+1
    self.model.maximize(sum((C-w)*self.x[i] for i,w in enumerate(df.get_column(weight))))
    self.solver=cp_model.CpSolver()
    self.status=self.solver.solve(self.model)
    if self.status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
      self.multiplets=PL.concat(df[i] for i in range(len(self.x)) if self.solver.value(self.x[i]))
      if verbose>=1:
        print(f'Maximum of objective function: {self.solver.objective_value}')
    else:
      self.multiplets=df.head(0)
      if verbose>=1:
        print("No solution found.")
    if verbose>=1:
      print(f"  status   : {self.solver.status_name(self.status)}")
      print(f"  conflicts: {self.solver.num_conflicts}")
      print(f"  branches : {self.solver.num_branches}")
      print(f"  wall time: {self.solver.wall_time} s")
    return self.multiplets
