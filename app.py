from dash import Dash, html, dash_table
import pandas as pd
import numpy as np
import visdcc
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px
from utils import return_hub_genes,return_degree_distribution, return_hub_genes_tf, return_connected_components, return_tf_projection, return_degree_distribution_graph_input, return_hub_tf, return_connected_components_projection
from dash.exceptions import PreventUpdate
import itertools



app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1'}]
)
clicks  = 1
clicks_2 = 1

column_names=["TF","Genes",'Edge','Score']
AD_GRN = pd.read_csv('../../../t2t_assembly/T2T_assembly_GRN/Stratification/COGDX/AD/AD_output_panda.txt', names = column_names,sep = ' ')
NCI_GRN = pd.read_csv('../../../t2t_assembly/T2T_assembly_GRN/Stratification/COGDX/NCI/output_panda_Control.txt', names = column_names,sep = ' ')


AD_GRN.sort_values(by='Score', ascending = False, inplace=True)
NCI_GRN.sort_values(by='Score', ascending = False, inplace=True)
AD_GRN = AD_GRN[AD_GRN.Score >= 9.00]
NCI_GRN = NCI_GRN[NCI_GRN.Score >=9.00]
AD_GRN.reset_index(drop=True, inplace= True)
edges_AD = []
NCI_GRN.reset_index(drop=True, inplace= True)
for row in AD_GRN.to_dict(orient = 'records'):
        source, target = row['TF'], row['Genes']
        edges_AD.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            'color': 'red',
        })
        break

node_list = list(set(AD_GRN['TF'].unique().tolist() + AD_GRN['Genes'].unique().tolist() + NCI_GRN['TF'].unique().tolist() + NCI_GRN['Genes'].unique().tolist() ))
#
TF_node_list = list(set(AD_GRN['TF'].unique().tolist()+ NCI_GRN['TF'].unique().tolist()))
nodes = [({'id': node_name, 'label': node_name, 'shape': 'dot', 'size':7}) 
         if node_name not in TF_node_list
         else
         ({'id': node_name, 'label': node_name, 'square': 'rec', 'size':5}) 
    
         for i,node_name in enumerate(node_list)]

'AD_TF_TF'

app.layout = html.Div([
    html.Div([html.H3('Select Network'),dcc.Dropdown(id = 'GRN_type',
                   options = [{'label': 'AD_GRN', 'value': 'AD_GRN'},
                              {'label': 'NCI_GRN', 'value': 'NCI_GRN'},
                              {'label': 'AD_TF_TF', 'value': 'AD_TF_TF'} ,
                              {'label': 'NCI_TF_TF', 'value': 'NCI_TF_TF'},
                              {'label': 'AD_Gene_Gene', 'value': 'AD_Gene_Gene'} ,
                              {'label': 'NCI_Gene_Gene', 'value': 'NCI_Gene_Gene'} 
                              ], value = 'AD_GRN')]),
    html.Br(),
    html.H3('Number of Edges'),  
    dcc.RangeSlider(0,4000,100, value = [0,100],id='edges-range-slider'), 

    
    visdcc.Network(id = 'net',
                   data = {'nodes': nodes, 'edges': edges_AD},
                   options = dict(height='600px', width = '80%') 
                   ),
    html.Br(),   
    html.Br(),
    dbc.Row(
            [
                dbc.Col(html.Button('Export Connected Components', id='export-all-component'), width={  "offset": 3}),
                dcc.Download(id="download-export-xlsx"),dbc.Col(html.Button('Export Edges', id='export-all-edges'), width={  "offset": 1}),
                dcc.Download(id="download-export-edges-xlsx")
            ]
        ),
    html.Br(),   
    html.Br(),
    dbc.Row(
            [
                dbc.Col(html.Div("Degree Distribution Genes"), width={  "offset": 2}),
                dbc.Col(html.Div("Degree Distribution TF"),width={  "offset": 2})
            ]
        ),
    dbc.Row(
            [
                dbc.Col(dcc.Graph(id = 'degree-genes', figure = {}, className = 'bar-chart')),
                dbc.Col(dcc.Graph(id = 'degree-tf', figure = {}, className = 'bar-chart')),
            ]
        ),

    html.Br(),   
    html.Br(),
    dbc.Row(
            [
                dbc.Col(html.Div("Hub Genes"), width={  "offset": 3}),
                dbc.Col(html.Div("Hub TF"), width={  "offset": 3})
            ]
        ),
    html.Br(),   
    html.Br(),
    dbc.Row(
            [
                dbc.Col(dash_table.DataTable(
                        id='hub-genes-paging',
                        columns=[
                        {"name": i, "id": i} for i in ['Index','Genes','Degree']
                        ],
                        page_current=0,
                        page_size=5,
                        page_action='custom'
                        ), width={ "size": 4, "offset": 1}),
                dbc.Col(dash_table.DataTable(
                        id='tf-genes-paging',
                        columns=[
                        {"name": i, "id": i} for i in ['Index','TF','Degree']
                        ],
                        page_current=0,
                        page_size=5,
                        page_action='custom'
                        ),width={ "size": 4, "offset": 2}),
            ]
        ),
    html.Div([html.H3('Select Nodes'),
        dcc.Dropdown(id = 'select_nodes',
            options = node_list,multi=True
                              , value = 'None')]),
    
    visdcc.Network(id = 'selected-net',
                   data = {'nodes': nodes, 'edges': edges_AD},
                   options = dict(height='700px', width = '80%'),style={"margin-left": "20px"} 
                   ),
      
])
#, color='#00ff00'
@app.callback(
        
    Output('net','options'),
    Output('net','data'),
    Output('selected-net','data'),
    Output('degree-genes','figure'),
    Output('degree-tf','figure'),
    Output('hub-genes-paging', 'data'),
    Output('tf-genes-paging', 'data'),
    Input('GRN_type','value'),
    Input('edges-range-slider', 'value'),
    Input('select_nodes', 'value'),
    Input('hub-genes-paging', "page_current"),
    Input('hub-genes-paging', "page_size"),
    Input('tf-genes-paging', "page_current"),
    Input('tf-genes-paging', "page_size"),
     )
def myfun(x, number_of_edges,selected_nodes, page_current_hub_genes,page_size_hub_genes, page_current_hub_tf,page_size_hub_tf):
    #print(number_of_edges)
    #print(edges_AD[number_of_edges[0]:number_of_edges[1]])
    #return {'nodes': nodes, 'edges': x}
    edges_AD = []
    edges_NCI = []
    edges_AD_selected = []
    edges_NCI_selected = []
    edges_selected = []
    node_list_selected = []

    if x == 'AD_GRN':


        for row in AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_dict(orient = 'records'):
            source, target = row['TF'], row['Genes']
            edges_AD.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            'color': 'red',
            })
        s = 0
        
        deg_genes, deg_tf = return_degree_distribution(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
        genes_degree_df, tf_degree_df = return_hub_genes_tf(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])

        for row in AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_dict(orient = 'records'):
            

            if len(selected_nodes)<1:
                break
            if row['TF'] in selected_nodes or row['Genes'] in selected_nodes:
                source, target = row['TF'], row['Genes']
                edges_AD_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                'color': 'red',
                })

                node_list_selected.append(row['TF'])
                node_list_selected.append(row['Genes']) 
            
        #print(edges_AD_selected)
        node_list = list(set(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],0].unique().tolist() + AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],1].unique().tolist())) 
        nodes = [({'id': node_name, 'label': node_name, 'shape': 'dot', 'size':7}) 
             if node_name not in TF_node_list
             else
             ({'id': node_name, 'label': node_name, 'square': 'rec', 'size':5})  
             for i,node_name in enumerate(node_list)]
        
        edges = edges_AD
        edges_selected = edges_AD_selected

    
    
    elif x == 'NCI_GRN':

        for row in NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_dict(orient = 'records'):
            source, target = row['TF'], row['Genes']
            edges_NCI.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            })


        deg_genes, deg_tf = return_degree_distribution(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
        genes_degree_df, tf_degree_df = return_hub_genes_tf(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
        

        for row in NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_dict(orient = 'records'):
            if len(selected_nodes)<1:
                break
            if row['TF'] in selected_nodes or row['Genes'] in selected_nodes:
                source, target = row['TF'], row['Genes']
                edges_NCI_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                })
                node_list_selected.append(row['TF'])
                node_list_selected.append(row['Genes']) 

        node_list = list(set(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],0].unique().tolist() + NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],1].unique().tolist())) 
        nodes = [({'id': node_name, 'label': node_name, 'shape': 'dot', 'size':7}) 
             if node_name not in TF_node_list
             else
             ({'id': node_name, 'label': node_name, 'square': 'rec', 'size':5})  
             for i,node_name in enumerate(node_list)]
        
        edges = edges_NCI
        edges_selected = edges_NCI_selected
    


    elif x == 'AD_TF_TF':
        TF_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
        edges_AD = []
        for row in TF_AD_projection.edges:
            source, target = row[0], row[1]
            edges_AD.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            }) 


        deg_tf = return_degree_distribution_graph_input(TF_AD_projection)
        genes_degree_df = pd.DataFrame(columns=['Index','Genes','Degree'])
        tf_degree_df = return_hub_tf(TF_AD_projection)
        print(selected_nodes)
        for row in TF_AD_projection.edges:
            if len(selected_nodes)<1:
                break
            if row[0] in selected_nodes or row[1] in selected_nodes:
                source, target = row[0], row[1]
                edges_AD_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                })
                node_list_selected.append(row[0])
                node_list_selected.append(row[1]) 

        
        node_list = list(set(list(itertools.chain(*list(TF_AD_projection.edges))))) 
        nodes = [({'id': node_name, 'label': node_name, 'size':5}) for i,node_name in enumerate(node_list)]
        edges = edges_AD
        edges_selected = edges_AD_selected
   
    elif x == 'AD_Gene_Gene':
        Gene_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:],  projection = 'Gene')
        edges_AD_selected = []
        edges_AD = []
        for row in Gene_AD_projection.edges:
            source, target = row[0], row[1]
            edges_AD.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            }) 


        deg_genes = return_degree_distribution_graph_input(Gene_AD_projection)
        genes_degree_df = return_hub_genes(Gene_AD_projection) 
        tf_degree_df = pd.DataFrame(columns=['Index','TF','Degree'])
        for row in Gene_AD_projection.edges:
            if len(selected_nodes)<1:
                break
            if row[0] in selected_nodes or row[1] in selected_nodes:
                source, target = row[0], row[1]
                edges_AD_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                })
                node_list_selected.append(row[0])
                node_list_selected.append(row[1]) 


        node_list = list(set(list(set(list(itertools.chain(*list(Gene_AD_projection.edges)))))))
        nodes = [({'id': node_name, 'label': node_name, 'size':5}) for i,node_name in enumerate(node_list)]
        edges = edges_AD
        edges_selected = edges_AD_selected


    elif x == 'NCI_TF_TF':
        TF_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
        edges_NCI = []
        for row in TF_NCI_projection.edges:
            source, target = row[0], row[1]
            edges_NCI.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            }) 
        

        deg_tf = return_degree_distribution_graph_input(TF_NCI_projection)
        genes_degree_df = pd.DataFrame(columns=['Index','Genes','Degree'])
        tf_degree_df = return_hub_tf(TF_NCI_projection)
        for row in TF_NCI_projection.edges:
            if len(selected_nodes)<1:
                break
            if row[0] in selected_nodes or row[1] in selected_nodes:
                source, target = row[0], row[1]
                edges_NCI_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                })
                node_list_selected.append(row[0])
                node_list_selected.append(row[1]) 


        node_list = list(set(list(set(list(itertools.chain(*list(TF_NCI_projection.edges))))) )) 
        nodes = [({'id': node_name, 'label': node_name, 'size':5}) for i,node_name in enumerate(node_list)]
        edges = edges_NCI
        edges_selected = edges_NCI_selected

    elif x == 'NCI_Gene_Gene':
        Gene_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:],  projection = 'Gene')
        edges_NCI_selected = []
        
        edges_AD = []
        for row in Gene_NCI_projection.edges:
            source, target = row[0], row[1]
            edges_NCI.append({
            'id': source +'__' + target,
            'from': source,
            'to': target,
            'width': 2,
            }) 


        deg_genes = return_degree_distribution_graph_input(Gene_NCI_projection)
        genes_degree_df = return_hub_genes(Gene_NCI_projection) 
        tf_degree_df = pd.DataFrame(columns=['Index','TF','Degree'])
        for row in Gene_NCI_projection.edges:
            if len(selected_nodes)<1:
                break
            if row[0] in selected_nodes or row[1] in selected_nodes:
                source, target = row[0], row[1]
                edges_NCI_selected.append({
                'id': source +'__' + target,
                'from': source,
                'to': target,
                'width': 2,
                })
                node_list_selected.append(row[0])
                node_list_selected.append(row[1]) 


        node_list = list(set(list(set(list(itertools.chain(*list(Gene_NCI_projection.edges)))))))
        nodes = [({'id': node_name, 'label': node_name, 'size':5}) for i,node_name in enumerate(node_list)]
        edges = edges_NCI
        edges_selected = edges_NCI_selected

    if 'GRN' in x:
        node_list_selected = list(set(node_list_selected))
        nodes_selected = [({'id': node_name, 'label': node_name, 'shape': 'dot', 'size':7}) 
             if node_name not in TF_node_list
             else
             ({'id': node_name, 'label': node_name, 'square': 'rec', 'size':5})  
             for i,node_name in enumerate(node_list_selected)]
    elif 'TF' in x:
        node_list_selected = list(set(node_list_selected))
        nodes_selected = [({'id': node_name, 'label': node_name, 'shape': 'rec', 'size':5}) for i,node_name in enumerate(node_list_selected)]
    
    elif 'Gene' in x:
        print('*'*100)
        node_list_selected = list(set(node_list_selected))
        nodes_selected = [({'id': node_name, 'label': node_name,  'shape': 'rec','size':5}) for i,node_name in enumerate(node_list_selected)]

    
    if 'GRN' in x:

        
        barchart_genes = px.bar(data_frame = deg_genes, labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
    
        barchart_tf = px.bar(data_frame = deg_tf, labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
        
        genes_degree_df= genes_degree_df.iloc[
        page_current_hub_genes*page_size_hub_genes:(page_current_hub_genes+ 1)*page_size_hub_genes
        ].to_dict('records')

        tf_degree_df= tf_degree_df.iloc[
        page_current_hub_tf*page_size_hub_tf:(page_current_hub_tf+ 1)*page_size_hub_tf
        ].to_dict('records')


    elif 'TF' in x:
        barchart_genes = px.bar(data_frame = pd.DataFrame(), labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
    
        barchart_tf = px.bar(data_frame = deg_tf, labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
        
        genes_degree_df= genes_degree_df.iloc[
        page_current_hub_genes*page_size_hub_genes:(page_current_hub_genes+ 1)*page_size_hub_genes
        ].to_dict('records')

        tf_degree_df= tf_degree_df.iloc[
        page_current_hub_tf*page_size_hub_tf:(page_current_hub_tf+ 1)*page_size_hub_tf
        ].to_dict('records')

    elif 'Gene' in x:
        barchart_tf = px.bar(data_frame = pd.DataFrame(), labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
    
        barchart_genes = px.bar(data_frame = deg_genes, labels={
                      "index":"Degree",
                     "value":"Number of Nodes"
                 },)
        
        genes_degree_df= genes_degree_df.iloc[
        page_current_hub_genes*page_size_hub_genes:(page_current_hub_genes+ 1)*page_size_hub_genes
        ].to_dict('records')

        tf_degree_df= tf_degree_df.iloc[
        page_current_hub_tf*page_size_hub_tf:(page_current_hub_tf+ 1)*page_size_hub_tf
        ].to_dict('records')
    
    #print(nodes_selected)
    
    
    
    
    return {'nodes': {'color': '#00ff00'}},{'nodes': nodes, 'edges': edges},{'nodes': nodes_selected, 'edges': edges_selected}, barchart_genes, barchart_tf,genes_degree_df,tf_degree_df





@app.callback(
    Output("download-export-xlsx", "data"),
    Input('GRN_type','value'),
    Input('edges-range-slider', 'value'),
    Input('export-all-component', 'n_clicks'),
    )
def clicks(x, number_of_edges,n_clicks):
    
    global clicks 
    
    if n_clicks is None or clicks == n_clicks:
        raise PreventUpdate

    else:
        
        if x == 'AD_GRN':
            name = 'AD_GRN_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            connected_components = return_connected_components(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv, name)
            
        elif x == 'NCI_GRN':
            name = 'NCI_GRN_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            connected_components = return_connected_components(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:])
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv, name)
        
        elif x == 'AD_TF_TF':
            name = 'AD_TF_TF_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            TF_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'TF')

            connected_components = return_connected_components_projection(TF_AD_projection)
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv, name)
        
        elif x == 'NCI_TF_TF':
            name = 'NCI_TF_TF_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            TF_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'TF')

            connected_components = return_connected_components_projection(TF_NCI_projection)
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv, name)
        
        elif x == 'AD_Gene_Gene':
            name = 'AD_Gene_Gene_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            Gene_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'Gene')

            connected_components = return_connected_components_projection(Gene_AD_projection)
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv,name)
        
        elif x == 'NCI_Gene_Gene':
            name = 'NCI_Gene_Gene_connected_components_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            Gene_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'Gene')

            connected_components = return_connected_components_projection(Gene_NCI_projection)
            clicks = n_clicks
            return dcc.send_data_frame(connected_components.to_csv, name)



@app.callback(
    Output("download-export-edges-xlsx", "data"),
    Input('GRN_type','value'),
    Input('edges-range-slider', 'value'),
    Input('export-all-edges', 'n_clicks'),
    )
def clicks_2(x, number_of_edges,n_clicks_2):
    
    global clicks_2 
    
    if n_clicks_2 is None or clicks_2 == n_clicks_2:
        raise PreventUpdate

    else:
        
        if x == 'AD_GRN':
            name = 'AD_GRN_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            clicks_2 = n_clicks_2
            return dcc.send_data_frame(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_csv, name)
            
        elif x == 'NCI_GRN':
            name = 'NCI_GRN_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            clicks_2 = n_clicks_2
            return dcc.send_data_frame(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:].to_csv, name)
        
        elif x == 'AD_TF_TF':
            TF_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'TF')

            x= []
            name = 'TF_AD_projection_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            for i in TF_AD_projection.edges:
                x.append([i[0],i[1]])
            return dcc.send_data_frame(pd.DataFrame(x).to_csv, name)  

        elif x == 'NCI_TF_TF':
            TF_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'TF')

            x= []
            name = 'TF_NCI_projection_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            for i in TF_NCI_projection.edges:
                x.append([i[0],i[1]])
            return dcc.send_data_frame(pd.DataFrame(x).to_csv, name)   

        elif x == 'AD_Gene_Gene':
            gene_AD_projection = return_tf_projection(AD_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'Gene')

            x= []
            name = 'gene_AD_projection_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'

            for i in gene_AD_projection.edges:
                x.append([i[0],i[1]])
            return dcc.send_data_frame(pd.DataFrame(x).to_csv, name)  

        elif x == 'NCI_Gene_Gene':
            gene_NCI_projection = return_tf_projection(NCI_GRN.iloc[number_of_edges[0]:number_of_edges[1],:], projection = 'Gene')

            x= []
            name = 'gene_NCI_projection_Edges_'+ str(number_of_edges[0])+'_' + str(number_of_edges[1])+'.csv'
            for i in gene_NCI_projection.edges:
                x.append([i[0],i[1]])
            return dcc.send_data_frame(pd.DataFrame(x).to_csv, name)       
        
    

if __name__ == '__main__':
    app.run(debug=True, port = 8020)