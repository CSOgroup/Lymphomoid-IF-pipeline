import os
import argparse
import numpy as np
import pandas as pd
import anndata
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")

## Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('--main_dir', type=str, required=True) # Absolute path to your main directory
parser.add_argument('--lymphomoids_to_process', type=str, nargs='+', default='all') # Can be 'all' or any space-separated list of lympomoid acquisitions under the folder 'Classified_cells_tables/', e.g. HLS09_s01_acq01_Ibrutinib01 HLS09_s01_acq01_Ibrutinib01
parser.add_argument('--celltype_query', type=str, default='Bcells') # Can be 'Bcells' (default), 'CD4', 'CD8' or 'Macrophages'
parser.add_argument('--celltype_query_proliferation_status', type=str, default='any') # Can be 'any' (default, all cells), 'proliferating', or 'not_proliferating'
parser.add_argument('--consider_othercells', type=str, default='False') # Whether or not to consider unassigned otherCells in the neighbourhood computation
parser.add_argument('--plot_delaunay', type=str, default='True') # Whether or not to plot the network
parser.add_argument('--custom_file_prefix', type=str, default='') # Custom file prefix for output (in case you run multiple times wih)
args = parser.parse_args()

this_palette = {"Bcells":"#00EE00","otherCell":"gray","Macrophages":"#EE00EE","CD4":"#FFD700","CD8":"#00CDCD"}

if args.lymphomoids_to_process[0]=='all':
    lymphomoids_to_process = os.listdir(args.main_dir+'Classified_cells_tables')
    lymphomoids_to_process = np.unique([x[6:len(x)-13] for x in lymphomoids_to_process if x.startswith("Table")])
else:
    lymphomoids_to_process = args.lymphomoids_to_process

if not os.path.exists(args.main_dir+'/Neighbourhood_analyses/'):
    os.makedirs(args.main_dir+'/Neighbourhood_analyses/')

consider_othercells = args.consider_othercells == 'True'
plot_delaunay = args.plot_delaunay == 'True'

## First, compute neighbourhood and distances acquisition by acquisition 
if consider_othercells:
    nei_df_all = pd.DataFrame(columns=[ "Bcells","CD4","CD8","Macrophages","otherCell","ll" ])
else:
    nei_df_all = pd.DataFrame(columns=[ "Bcells","CD4","CD8","Macrophages","ll" ])
dist_df_all = nei_df_all.copy()
df_all = pd.DataFrame(columns=[ "CellType_marker","ll" ])

for ll in lymphomoids_to_process:
    print("\n")
    print("------ Processing ", ll," ------")
    # Reading and preprocessing table
    df = pd.read_table( args.main_dir+'Classified_cells_tables/Table_'+ll+'_AllCells.txt' )
    # df = pd.read_table( '/mnt/ndata/daniele/elisa_lymphomoids/Processed/Pipeline_test/Classified_cells_tables/Table_'+ll+'_AllCells.txt' )
    
    if not consider_othercells:
        df = df.loc[df.CellType_marker != "otherCell",:]
    df.CellType_marker = df.CellType_marker.str.replace(" ","")
    df.CellType_marker = df.CellType_marker.str.replace("Tcells","")
    if (np.sum(df.CellType_marker==args.celltype_query)<5):
        print("Skipping", ll,", less than 5 query",args.celltype_query)
        continue
    adata = anndata.AnnData(df.iloc[:,3:7],obsm={"spatial": np.array(df.iloc[:,7:9])},dtype=np.float32)
    adata.obs['CellType_color'] = df.CellType_color.values
    adata.obs['CellType_marker'] = df.CellType_marker.values
    adata.obs['is_proliferating'] = df.is_proliferating.values
    sq.gr.spatial_neighbors(adata,coord_type='generic',delaunay=True)
    plt.style.use('dark_background')
    if plot_delaunay:
        sc.pl.spatial(adata, neighbors_key="spatial_neighbors", spot_size=5, color='CellType_marker', palette=this_palette, edges=True, edges_width=0.05, edges_color='#EBEBEB', show=False)
        plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_delaunay_plot_beforeFiltering.pdf')
    conns, dists = adata.obsp["spatial_connectivities"], adata.obsp["spatial_distances"]
    distance_percentile = 95.0
    threshold = np.percentile(np.array(dists[dists != 0]).squeeze(), distance_percentile)
    conns[dists > threshold] = 0
    dists[dists > threshold] = 0
    conns.eliminate_zeros()
    dists.eliminate_zeros()
    adata.obs['is_boundary'] = "no"
    cells_boundary = adata.obs.loc[adata.obsp['spatial_connectivities'].sum(axis=1)<4,].index
    adata.obs.loc[cells_boundary,'is_boundary'] = "yes"
    if plot_delaunay:
        sc.pl.spatial(adata, neighbors_key="spatial_neighbors", spot_size=5, color='CellType_marker', palette=this_palette, edges=True, edges_width=0.05, edges_color='#EBEBEB', show=False)
        plt.tight_layout()
        plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_delaunay_plot_afterFiltering.pdf')
        sc.pl.spatial(adata, neighbors_key="spatial_neighbors", spot_size=5, color='is_boundary', palette={"yes":"red","no":"#EBEBEB"}, edges=True, edges_width=0.05, edges_color='#EBEBEB', show=False)
        plt.tight_layout()
        plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_delaunay_plot_Boundary.pdf')
    if args.celltype_query_proliferation_status=="proliferating":
        query_cells = adata.obs.loc[(adata.obs['CellType_marker']==args.celltype_query) & (adata.obs['is_proliferating']),].index
        title_neighComp = ll+",\nneighbourhood of proliferating "+args.celltype_query
    elif (args.celltype_query_proliferation_status=="not_proliferating"):
        query_cells = adata.obs.loc[(adata.obs['CellType_marker']==args.celltype_query) & (~adata.obs['is_proliferating']),].index
        title_neighComp = ll+",\nneighbourhood of non proliferating "+args.celltype_query
    else:
        query_cells = adata.obs.loc[(adata.obs['CellType_marker']==args.celltype_query),].index
        title_neighComp = ll+",\nneighbourhood of "+args.celltype_query
    query_cells = pd.Index(list(set(query_cells)-set(cells_boundary)))
    if (len(query_cells)<5):
        print("Skipping", ll,", less than 5 query",args.celltype_query,"after subsetting for proliferation status")
        continue
    if consider_othercells:
        nei_df = pd.DataFrame(np.NaN,query_cells,[ "Bcells","CD4","CD8","Macrophages","otherCell" ])
        this_df = df.loc[:,["CellType_marker","is_proliferating"] ]
        this_df['ll'] = ll
        this_df = this_df.drop(columns="is_proliferating")
        df_all = pd.concat([df_all,this_df],axis=0)
    else:
        nei_df = pd.DataFrame(np.NaN,query_cells,[ "Bcells","CD4","CD8","Macrophages" ])
        this_df = df.loc[df.CellType_marker!="otherCell",["CellType_marker","is_proliferating"] ]
        this_df['ll'] = ll
        this_df = this_df.drop(columns="is_proliferating")
        df_all = pd.concat([df_all,this_df],axis=0)
    dist_df = nei_df.copy()

    query_spatial_conn = adata.obsp['spatial_connectivities'][:,:].toarray()
    query_spatial_dist = adata.obsp['spatial_distances'][:,:].toarray()
    for q in tqdm(query_cells):
        nei = adata.obs.loc[query_spatial_conn[:,adata.obs_names==q]==1,].copy()
        dist = query_spatial_dist[ adata.obs_names.isin(nei.index) ,adata.obs_names==q]
        nei['distance'] = dist
        nei = nei.loc[:,[ "CellType_marker","distance" ]]
        tabb = nei.CellType_marker.value_counts()/nei.shape[0]
        nei_df.loc[q,tabb.index.tolist()] = tabb.tolist()
        tabb = nei.groupby(['CellType_marker']).mean()
        dist_df.loc[q,tabb.index.tolist()] = tabb['distance']

    # Plotting neigh composition
    nei_df_melt = pd.melt(nei_df)
    nei_df_melt.columns = [ "Neighbouring cell type","Relative frequency" ]
    plt.style.use('classic')
    plt.figure()
    sns.barplot(data=nei_df_melt, x="Neighbouring cell type", y="Relative frequency", palette=this_palette)
    plt.ylim(0, 1)
    plt.title( title_neighComp )
    plt.tight_layout()
    plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_NeighbourhoodFrequencies.pdf')
    plt.close()
    nei_df.to_csv(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_NeighbourhoodFrequencies_FullTable.csv')

    # Plotting distance density distributions
    dist_df_melt = pd.melt(dist_df)
    dist_df_melt = dist_df_melt.dropna()
    dist_df_melt.columns = [ "Neighbouring cell type","Distance (µm)" ]
    plt.figure()
    sns.kdeplot(data=dist_df_melt, hue="Neighbouring cell type", x="Distance (µm)",linewidth=2, palette=this_palette)
    plt.tight_layout()
    plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_DistanceDistributions.pdf')
    plt.close()
    dist_df.to_csv(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+ll+'_DistanceDistributions_FullTable.csv')

    # Concatenating into full dataframes
    nei_df['ll'] = ll
    dist_df['ll'] = ll
    nei_df_all = pd.concat([nei_df_all,nei_df],axis=0)
    dist_df_all = pd.concat([dist_df_all,dist_df],axis=0)

## Second, collapse lymphomoids
nei_df_all['PatientLymphomoidName'] = [thisll[::-1].partition('_')[0][::-1] for thisll in nei_df_all['ll']]
dist_df_all['PatientLymphomoidName'] = [thisll[::-1].partition('_')[0][::-1] for thisll in dist_df_all['ll']]
df_all['PatientLymphomoidName'] = [thisll[::-1].partition('_')[0][::-1] for thisll in df_all['ll']]
nei_df_all.to_csv(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'AllLymphomoids_NeighbourhoodFrequencies_FullTable.csv')
dist_df_all.to_csv(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'AllLymphomoids_DistanceDistributions_FullTable.csv')

# nei_df_all = pd.read_csv( "/mnt/ndata/daniele/elisa_lymphomoids/Processed/Pipeline_test/Neighbourhood_analyses/LP08Any_AllLymphomoids_NeighbourhoodFrequencies_FullTable.csv" )
nei_df_all.loc[nei_df_all.Bcells.isna(),"Bcells"] = 0
nei_df_all.loc[nei_df_all.CD4.isna(),"CD4"] = 0
nei_df_all.loc[nei_df_all.CD8.isna(),"CD8"] = 0
nei_df_all.loc[nei_df_all.Macrophages.isna(),"Macrophages"] = 0
# nei_df_all = nei_df_all.loc[nei_df_all.ll.isin(df_all.ll.unique()),]

# Plotting neigh composition
nei_df_all = nei_df_all.drop(columns="ll")
df_all = df_all.drop( columns = 'll' )
# nei_df_all = nei_df_all.drop(columns="Unnamed: 0")
dist_df_all = dist_df_all.drop(columns="ll")

nei_df_grouped = nei_df_all.groupby(['PatientLymphomoidName']).mean()
# nei_df_grouped.sum(1)

# Obs vs Exp
df_all_grouped = nei_df_grouped.copy()
for i in df_all_grouped.index:
    for j in df_all_grouped.columns:
        df_all_grouped.loc[i,j] = np.sum((df_all.CellType_marker==j) & (df_all.PatientLymphomoidName==i) )
df_all_grouped = df_all_grouped.div(df_all_grouped.sum(1),axis=0)
obs_vs_exp_log_ratio = np.log2(nei_df_grouped/df_all_grouped)
obs_vs_exp_log_ratio['PatientLymphomoidName'] = obs_vs_exp_log_ratio.index
oe = pd.melt(obs_vs_exp_log_ratio,'PatientLymphomoidName')
oe.columns = ["PatientLymphomoidName", "CellType","Obs_vs_Exp_log2_ratio" ]
plt.style.use('classic')
plt.figure()
sns.barplot(data=oe, x="PatientLymphomoidName", y="Obs_vs_Exp_log2_ratio", hue="CellType", palette=this_palette)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'Obs_vs_Exp_log2_ratio.pdf')
plt.close()
oe.to_csv(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'AllLymphomoids_Obs_vs_Exp_log2_ratio_Table.csv')

# Barplot overall
nei_df_grouped['PatientLymphomoidName'] = nei_df_grouped.index.values
plt.style.use('classic')
plt.figure()
fig, ax = plt.subplots()
groups = nei_df_grouped['PatientLymphomoidName'].values.tolist()
values_Bcells = nei_df_grouped['Bcells'].values.tolist()
values_CD4 = nei_df_grouped['CD4'].values.tolist()
values_CD8 = nei_df_grouped['CD8'].values.tolist()
values_Macrophages = nei_df_grouped['Macrophages'].values.tolist()
ax.bar(groups, values_Bcells, color = "#00EE00", edgecolor = "black", linewidth = 1,label='Bcells')
ax.bar(groups, values_CD4, bottom = values_Bcells, color = "#FFD700", edgecolor = "black", linewidth = 1,label='CD4')
ax.bar(groups, values_CD8, bottom = np.add(values_Bcells,values_CD4), color = "#00CDCD", edgecolor = "black", linewidth = 1,label='CD8')
ax.bar(groups, values_Macrophages, bottom = np.add(np.add(values_Bcells,values_CD4),values_CD8), color = "#EE00EE", edgecolor = "black", linewidth = 1,label='Macrophages')
if consider_othercells:
    values_otherCell = nei_df_grouped['otherCell'].values.tolist()
    ax.bar(groups, values_otherCell, bottom = np.add(np.add(np.add(values_Bcells,values_CD4),values_CD8),values_Macrophages), color = "gray", edgecolor = "black", linewidth = 1,label='otherCell')
ax.legend(loc='lower right')
plt.xticks(rotation=45)
plt.ylabel('Relative frequency' )
plt.tight_layout()
# plt.savefig("/mnt/ndata/daniele/elisa_lymphomoids/Processed/Pipeline_test/Neighbourhood_analyses/LP08Any_AllLymphomoids_NeighbourhoodFrequencies_temp.pdf" )
plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'AllLymphomoids_NeighbourhoodFrequencies.pdf')
plt.close()

for col in list(set(['Bcells','CD4','CD8','Macrophages','otherCell']) & set(dist_df_all.columns.values.tolist())):
    this_dist_df_all = dist_df_all.loc[:,[col,"PatientLymphomoidName"]].copy().dropna()
    plt.figure()
    sns.kdeplot(data=this_dist_df_all, hue='PatientLymphomoidName', x=col,linewidth=2,palette="Set2")
    plt.tight_layout()
    plt.savefig(args.main_dir+'/Neighbourhood_analyses/'+args.custom_file_prefix+'AllLymphomoids_'+col+'_DistanceDistributions.pdf')
    plt.close() 













