#Function foldchange_log2() to compute log2 foldchange of x and y

foldchange_log2<-function(x,y,logged2)
{
	if(!logged2){x=log2(x)}
	m1 <- rowMeans(x[,which(y==1)])
	m2 <- rowMeans(x[,which(y==2)])
	return(m2-m1)
}

#Function labeledge() to have a list of label for each edge
#it takes the name of database:"kegg","reactome","nci"
#and the name of pathway and gives a dataframe with
#the characteristics about the type of edge between genes

labeledge<-function(database,name)
{
	if(database=="kegg")
	{
		# p<-kegg[[name]]@edges
		p<-graphite::edges(kegg[[name]])
		src<-vector()
		dest<-vector()
		srcdest<-vector()
		for(i in 1:dim(p)[[1]])
		{
		  src[i]<-p[i,2]
		  dest[i]<-p[i,4]
			if(as.character(p[i,5])=="undirected")
			  srcdest[i]<-paste(src[i],dest[i],sep="<->") else
			  srcdest[i]<-paste(src[i],dest[i],sep="->")
		}
		# db1<-cbind(src,dest,srcdest,as.character(p[,3]),as.character(p[,4]))
		db1<-cbind(src,dest,srcdest,as.character(p[,5]),as.character(p[,6]))
		colnames(db1)<-c("src","dest","srcdest","link","type_edge")


	}
	if(database=="reactome")
	{
		p<-reactome[[name]]@edges
		src<-vector()
		dest<-vector()
		srcdest<-vector()

		for(i in 1:dim(p)[[1]])
		{
			src[i]<-strsplit(p[i,1],"EntrezGene:")[[1]][2]
			dest[i]<-strsplit(p[i,2],"EntrezGene:")[[1]][2]
			if(as.character(p[i,3])=="undirected")
				srcdest[i]<-paste(src[i],dest[i],sep="<->")
			else
				srcdest[i]<-paste(src[i],dest[i],sep="->")
		}
		db1<-cbind(src,dest,srcdest,as.character(p[,3]),as.character(p[,4]))
		colnames(db1)<-c("src","dest","srcdest","link","type_edge")

	}
	if(database=="nci")
	{
		p<-reactome[[name]]@edges
		src<-vector()
		dest<-vector()
		srcdest<-vector()
		for(i in 1:dim(p)[[1]])
		{
			src[i]<-strsplit(p[i,1],"EntrezGene:")[[1]][2]
			dest[i]<-strsplit(p[i,2],"EntrezGene:")[[1]][2]
			if(as.character(p[i,3])=="undirected")
				srcdest[i]<-paste(src[i],dest[i],sep="<->")
			else
				srcdest[i]<-paste(src[i],dest[i],sep="->")
		}
		db1<-cbind(src,dest,srcdest,as.character(p[,3]),as.character(p[,4]))
		colnames(db1)<-c("src","dest","srcdest","link","type_edge")
	}
	db1
}

#Function src_dest() that gives a dataframe of source nodes and destination node from 
# a adjacency matrix

src_dest<-function(adjmatr)
{
	src<-vector()
	dest<-vector()
	srcdest<-vector()
	c<-1
	for(i in 1:dim(adjmatr)[[1]])
	{
		for(j in 1:dim(adjmatr)[[2]])
		{
			if(adjmatr[i,j]==1 && adjmatr[j,i]==0)
			{
				src[c]<-rownames(adjmatr)[i]
				dest[c]<-colnames(adjmatr)[j]	
				srcdest[c]<-paste(src[c],dest[c],sep="->")
				c<-c+1
			}
			else if(adjmatr[i,j]==1 && adjmatr[j,i]==1)
			{
				src[c]<-rownames(adjmatr)[i]
				dest[c]<-colnames(adjmatr)[j]	
				srcdest[c]<-paste(src[c],dest[c],sep="<->")
				c<-c+1
			}

		}
	}
	db<-cbind(src,dest,srcdest)
	colnames(db)<-cbind("scr","dest","scrdest")
	db
}

#Function find_positions(), it accepts two vector and find the position of elements of b in a

find_positions<-function(a,b)
{
	pos<-NULL
	for(i in 1:length(b))
	{
		pos<-c(pos,which(a==b[i]))
	}
	pos
}


#Function merge_pathway() to merge graphs
#The function requires the name of pathways and the relative database.
#The function transforms the pathway in graph and the performs the merging of the relative
#adjacency matrices
#the function returns the adjacency matrix after all merging.

merge_pathway<-function(db,name_path)
{
	if(db=="kegg")
	{
		p_entrez_kegg<-list()
		for(i in 1:length(name_path))
		{
			p_entrez_kegg[[i]]<-convertIdentifiers(kegg[[name_path[[i]]]],"entrez")
		
		}
		p_entrez_gr<-list()
		for(i in 1:length(p_entrez_kegg))
		{
			p_entrez_gr[[i]]<-pathwayGraph(p_entrez_kegg[[i]])
		}
		#Union matrix
		union_matr<-as(p_entrez_gr[[1]],"matrix")
		if(length(name_path)>1)
		{
			for(i in 2:(length(name_path)))
			{
				union_matr<-merge_adj_matr(union_matr,as(p_entrez_gr[[i]],"matrix"),0)
			}
		}
	}
	if(db=="reactome")
	{
		p_entrez_reac<-list()
		for(i in 1:length(name_path))
		{
			p_entrez_reac[[i]]<-convertIdentifiers(reactome[[name_path[[i]]]],"entrez")
		
		}
		p_entrez_gr<-list()
		for(i in 1:length(p_entrez_reac))
		{
			p_entrez_gr[[i]]<-pathwayGraph(p_entrez_reac[[i]])
		}
	
		#Union matrix
		union_matr<-as(p_entrez_gr[[1]],"matrix")
		if(length(name_path)>1)
		{
			for(i in 2:(length(name_path)))
			{
				union_matr<-merge_adj_matr(union_matr,as(p_entrez_gr[[i]],"matrix"),0)
			}
		}
	}
	if(db=="nci")
	{
		p_entrez_nci<-list()
		for(i in 1:length(name_path))
		{
			p_entrez_nci[[i]]<-convertIdentifiers(nci[[name_path[[i]]]],"entrez")
		
		}
		p_entrez_gr<-list()
		for(i in 1:length(p_entrez_kegg))
		{
			p_entrez_gr[[i]]<-pathwayGraph(p_entrez_nci[[i]])
		}
	
		#Union matrix
		union_matr<-as(p_entrez_gr[[1]],"matrix")
		if(length(name_path)>1)
		{
			for(i in 2:(length(name_path)))
			{
				union_matr<-merge_adj_matr(union_matr,as(p_entrez_gr[[i]],"matrix"),0)
			}
		}
	}
	union_matr
}

#Function graph_DE(), it calculates all shortest paths between pairs of vertices.
#More precisely, between the from vertex to the vertices given in to. It uses a breadth-first
#search for unweighted graphs and Dijkstra?s algorithm for weighted ones.
#The latter only supports non-negative edge weights.
#So I wish to find a path for each gene present in the list of DE.
#I need to perform two tytpe of operations: 1) find the paths for each DEG; 2)
#transform the pos in name.

#graph_path is the union of all significant pathways found (graphNEL object)
#graph_micro must be a graphNEL object
#id_DE is the identificative of DE found

#the function returns: 
#1) a "igraph" object relative to the connection to DE through microarray genes
#2) vertex names of the previous graph
#3) a list with the shortest paths between the De
#4) type of edges

graph_DE<-function(graph_path,graph_micro,id_DE,type)
{
	#Here I find the position of DEGs on the graph relative
	#to pathway containing only microarray genes
	pos_sig_graph<-find_positions(graph_micro@nodes,id_DE)
	sg_igraph<-igraph::igraph.from.graphNEL(graph_micro)
	#here I find the shortesh paths that start from a DEGs and arrive to a DEGs
	gene_path<-list()
	for(i in 1:length(pos_sig_graph))
	{
		gene_path[[i]]<-igraph::get.shortest.paths(sg_igraph,pos_sig_graph[i],pos_sig_graph,mode=type,output="vpath")
		gene_path[[1]]$vpath[[2]]
	}
	#One of the possible solution to sinthesize data is to obtain a list 
	#of all genes of interest and then to extract information of only those genes
	#Here I take all genes found in the shortest path computed on the subgraph 
	#of the pathway containing only microarray genes
	#Then I extract the subgraph of that genes
	all_interested_genes<-NULL
	for(i in 1:length(gene_path))
	{
		for(j in 1:length(gene_path[[i]]$vpath))
		{
			all_interested_genes<-c(all_interested_genes,as_ids(gene_path[[i]]$vpath[[j]]))
		}
	}
	all_interested_genes<-unique(all_interested_genes)
	sg_de<-subGraph(intersect(all_interested_genes,graph_path@nodes),graph_path)
	sg_igraph_de<-igraph::igraph.from.graphNEL(sg_de)
	path<-list()
	c<-1
	for(i in 1:length(gene_path))
	{
		for(j in 1:length(gene_path[[i]]$vpath))
		{
			if(length(gene_path[[i]]$vpath[[j]])!=0)
			{
				path[[c]]<-gene_path[[i]]$vpath[[j]]
				c<-c+1
			}
		}
	}
	#This will be the list of nodes
	#Then I need the list of edges path
	output<-list(sg_igraph_de,sg_de@nodes,path)
	output
}

#Function assignlabel() that given a scrdest find the type of link for links indicated 
#Arguments
#The database with information about the edges->dbtot
#the scrdest on which you want to obtain the description

assignlabel<-function(db_tot,db)
{
	label<-vector()
	for(i in 1:dim(db)[[1]])
	{
		pos<-NULL
		pos<-which(db_tot[,3]==db[i,3])
		pos1<-NULL
		pos1<-which(db_tot[,3]==paste(strsplit(db[i,3],'<->')[[1]][2],strsplit(db[i,3],'<->')[[1]][1],sep='<->'))
		pos2<-NULL
		pos2<-which(db_tot[,3]==paste(strsplit(db[i,3],'<->')[[1]][2],strsplit(db[i,3],'<->')[[1]][1],sep='->'))
		pos3<-NULL
		pos3<-which(db_tot[,3]==paste(strsplit(db[i,3],'<->')[[1]][1],strsplit(db[i,3],'<->')[[1]][2],sep='->'))
		if(length(pos)==0 && length(pos1)==0 && length(pos2)==0 && length(pos3)==0)
		{
			label[i]<-c("coexpression")
		}
		else if(length(pos)==0 && length(pos1)==0 && length(pos2)==0 && length(pos3)!=0)
		{
			label[i]<-db_tot[pos3,5]
		}
		else if(length(pos)==0 && length(pos1)==0 && length(pos2)!=0 && length(pos3)==0)
		{
			label[i]<-db_tot[pos2,5]
		}
		else if(length(pos)==0 && length(pos1)==0 && length(pos2)!=0 && length(pos3)!=0)
		{
			label[i]<-db_tot[pos2,5]
			db[i,3]<-paste(strsplit(db[i,3],'<->')[[1]][1],strsplit(db[i,3],'<->')[[1]][2],sep='-><-')	
		}
		else if(length(pos)==0 && length(pos1)!=0 && length(pos2)==0 && length(pos3)==0)
		{
			label[i]<-db_tot[pos1,5]
		}
		else if(length(pos)==0 && length(pos1)!=0 && length(pos2)!=0 && length(pos3)==0)
		{
			label[i]<-paste(db_tot[pos2,5],",",db_tot[pos1,5])
			db[i,3]<-paste(strsplit(db[i,3],'<->')[[1]][1],strsplit(db[i,3],'<->')[[1]][2],sep='-><->')	
		}


		else if(length(pos)>1)
		{
			l<-db_tot[pos[1],5]
			for(m in 2:length(pos))
			{
				l<-paste(l,db_tot[pos[m],5],sep=",")
			}
		label[i]<-l
		}
		else
			label[i]<-db_tot[pos,5]
	}
	db_annot<-cbind(db,label)
	db_annot
}

#Function edge_paths() that gives back the list of edges starting from the list of nodes

edge_paths<-function(db_tot,path)
{
	edgepath<-list()
	app<-NULL
	for(i in 1:length(path))
	{
		for(j in 1:(length(path[[i]])-1))
		{
			app<-c(app,db_tot[intersect(which(db_tot[,1]==path[[i]][[j]]),which(db_tot[,2]==path[[i]][[j+1]])),3])
		}
		edgepath[[i]]<-app
		app<-NULL
	}
	edgepath
}

#The function inters_genes() returns a list on intersections  
#of each pathway with microarray not DE genes and intersections  
#of each pathway with DE

inters_genes<-function(name_path,name_micro,name_diff_genes,database)
{
	num_path<-length(name_path)
	inters_no_de<-list()
	inters_de<-list()
	for(i in 1:num_path)
	{
		inters_no_de[[i]]<-intersect(name_micro,rownames(merge_pathway(database,name_path[i])))
		inters_de[[i]]<-intersect(name_diff_genes,rownames(merge_pathway(database,name_path[i])))
	}
	names(inters_no_de)<-name_path
	names(inters_de)<-name_path
	
	inters<-list(inters_no_de,inters_de)
	inters
}

#Function new_path() that gives paths with  PIRSF

new_path<-function(DEGs,SUPF,list_nodes,list_edges)
{
	new_nodes<-NULL
	new_edges<-NULL
	new_edges_supf<-NULL
	j<-length(list_nodes)
	#I create a list of genes belonging to PIR superfamily
	list_SUPF<-NULL
	for(i in 1:length(SUPF))
	{
		list_SUPF<-c(list_SUPF,SUPF[[i]])
	}
	if(length(list_nodes)==2)
	{
		new_nodes<-list_nodes
		new_edges<-list_edges
		new_edges_supf<-list_edges
		return(list(new_nodes, new_edges,new_edges_supf))
	}
	else
	{
		new_nodes[1]<-list_nodes[1]
		c<-2
		for (d in 2:(j-1))
		{
			if (is.element(list_nodes[d],DEGs) && !is.element(list_nodes[d],new_nodes) && !is.element(list_edges[d-1],new_edges)) #We verify if the next nodes of X is a DEG
 			{
				new_nodes[c] = list_nodes[d]
				new_edges[c-1]= list_edges[d-1]
				c=c+1
			}
			if(is.element(list_nodes[d],list_SUPF))
			{
				for(s in 1:length(SUPF))
				{
					if (is.element(list_nodes[d],SUPF[[s]]))
					{
						new_node= names(SUPF)[s]
						new_edge= list_edges[d-1]
						if(!is.element(new_edge,new_edges))
						{
							new_nodes[c]<-new_node
							new_edges[c-1]<-new_edge	
							c=c+1
						}
					}
				}
			}
			else#The next nodes of X will be a not DEGs
			{			
				new_nodes[c] = list_nodes[d]
				new_edges[c-1] = list_edges[d-1]
				c=c+1
			}
		}
		new_nodes[j]<-list_nodes[j]
		new_edges[j-1]<-list_edges[j-1]
		c<-1
		j<-1
		i<-j+1
		while(c<=length(list_edges))
		{
			new_edges_supf[c]<-new_edges[c]
			while(j<=i)
			{
				new_edges_supf[c]<-sub(list_nodes[j],new_nodes[j],new_edges_supf[c])
				j<-j+1
			}
			i<-j
			j<-j-1
			c<-c+1
		}
		if(length(which(duplicated(new_nodes)==TRUE))!=0)
			new_nodes<-new_nodes[-which(duplicated(new_nodes)==TRUE)]
		if(length(which(duplicated(new_edges_supf)==TRUE))!=0)
		{
			new_edges<-new_edges[-which(duplicated(new_edges_supf)==TRUE)]
			new_edges_supf<-new_edges_supf[-which(duplicated(new_edges_supf)==TRUE)]
		}
		#Now to delete the connection in which there is the connection with himself
		pos<-NULL
		for(i in 1:length(new_nodes))
		{
			for(j in 1:length(new_edges_supf))
			{
				if(length(gregexpr(new_nodes[i],new_edges_supf[j],perl=T)[[1]])==2)
					pos<-c(pos,j)
			}
		}
		if(length(pos)!=0)
		{
			new_edges<-new_edges[-pos]
			new_edges_supf<-new_edges_supf[-pos]
		}
		return (list(new_nodes, new_edges,new_edges_supf))
	}
}

#Function plot.paths() that designs graph from the paths found
#plot.paths(DEG_focal[[1]],names(pirsf2),path_node_edge2)
#DEGs<-DEG_focal[[1]]
#SUPF<-names(pirsf2)
#paths<-path_node_edge2

plot.paths<-function(DEGs,SUPF,paths)
{
	#The first step is to merge nodes and edges
	#Merging of nodes
	nodes<-vector()
	for(i in 1:length(paths))
	{
		for(j in 1:length(paths[[i]][[1]]))
			nodes<-c(nodes,paths[[i]][[1]][j])
	}
	nodes<-unique(nodes)#list of nodes
	#Merging of edges
	edges<-vector()
	for(i in 1:length(paths))
	{
		for(j in 1:length(paths[[i]][[3]]))
			edges<-c(edges,paths[[i]][[3]][j])
	}
	edges<-unique(edges)#list of nodes
	#The best solution is to obtain an adjacency matrix
	#Step 1: Creation of the adjacency matrix
	adj_model<-matrix(rep(0,length(nodes)*length(nodes)),nrow=length(nodes),ncol=length(nodes))
	rownames(adj_model)<-colnames(adj_model)<-nodes
	#Step2 fill the adjacency matrix
	for(i in 1:length(rownames(adj_model)))
	{
		pos<-NULL
		for(j in 1:length(edges))
		{
			#I verify in which edge the node is present as source
			sep<-strsplit(edges[j],nodes[i],fixed=TRUE)
			if(length(sep[[1]])>1)
			{	
				#In which edge is present the nodes i as source
				pos<-c(pos,j)
				#I keep the info relative to which type of link is
				#and the node that receive the node
			}
		}
		if(length(pos)!=0)
		{
			for(j in 1:length(pos))
			{
				#Possible types of edges
				type1<-strsplit(edges[pos[j]],paste(nodes[i],"->",sep=""))
				type2<-strsplit(edges[pos[j]],paste(nodes[i],"<->",sep=""))
				type3<-strsplit(edges[pos[j]],paste(nodes[i],"-><-",sep=""))
				#Edge actually present
				type_present<-which(c(length(type1[[1]]),length(type2[[1]]),length(type3[[1]]))>1)
				if(length(type_present)==1)
				{
					if(type_present==1)
						adj_model[i,which(colnames(adj_model)==type1[[1]][2])]=1
					if(type_present==2)
					{
						adj_model[i,which(colnames(adj_model)==type2[[1]][2])]=100
						adj_model[which(colnames(adj_model)==type2[[1]][2]),i]=100
					}
				}
				if(length(type_present)>1)
				{
					adj_model[i,which(colnames(adj_model)==type3[[1]][2])]=1
					adj_model[which(colnames(adj_model)==type3[[1]][2]),i]=1
				}
			}
		}
	}
	#Now I have the adjacency matrix
	#The next step is to create the graph
	#Transform the adjacency matrix in igraph object
	#graphmodel<-graph.adjacency(adj_model)
	#vcol<-rep("yellow",length(V(graphmodel)$name))
	#vcol[find_positions(V(graphmodel)$name,DEGs)]<-"green"
	#vcol[find_positions(V(graphmodel)$name,SUPF)]<-"light blue"
	vcol<-rep("yellow",length(rownames(adj_model)))
	vcol[find_positions(rownames(adj_model),DEGs)]<-"green"
	vcol[find_positions(rownames(adj_model),SUPF)]<-"light blue"
	plotGraph(adj_model,vc=vcol,tcltk=FALSE,directed=TRUE)
	#plot(graphmodel,vertex.label=V(graphmodel)$name,vertex.color=vcol,vertex.size=20,vertex.label.cex=0.8)
	return(adj_model)
}

#Function block_matrix() creates a block matrix only for two matrices

block_matrix<-function(matr1,matr2)
{
	nr<-dim(matr1)[[1]]+dim(matr2)[[1]]
	nc<-dim(matr1)[[2]]+dim(matr2)[[2]]
	tot_matr<-matr<-matrix(c(rep(0,nr*nc)),nrow=nr,ncol=nc,byrow=T)
	for(i in 1:dim(matr1)[[1]])
	{
		for(j in 1:dim(matr1)[[2]])
		{
			tot_matr[i,j]=matr1[i,j]
		}
	}
	for(i in 1:dim(matr2)[[1]])
	{
		for(j in 1:dim(matr2)[[2]])
		{
			tot_matr[i+dim(matr1)[[1]],j+dim(matr1)[[2]]]=matr2[i,j]
		}
	}
	rownames(tot_matr)<-c(rownames(matr1),rownames(matr2))
	colnames(tot_matr)<-c(colnames(matr1),colnames(matr2))
	tot_matr
}

#db_tot is constitued by scr, dest, srcdest and label
#ls<-list of groups
assignlabelgroup<-function(db_tot,ls)
{
	srcg<-vector()
	destg<-vector()
	srcdestg<-vector()
	labelg<-vector()
	for(i in 1:dim(db_tot)[[1]])
	{
		app<-NULL
		for(j in 1:length(ls))
		{
			if(is.element(db_tot[i,1],ls[[j]]))
				app<-names(ls)[j]
		}
		if(length(app)==1)
			srcg[i]=app
		else
			srcg[i]=db_tot[i,1]
		labelg[i]<-db_tot[i,4]
	}
	for(i in 1:dim(db_tot)[[1]])
	{
		app<-NULL
		for(j in 1:length(ls))
		{
			if(is.element(db_tot[i,2],ls[[j]]))
				app<-names(ls)[j]
		}
		if(length(app)==1)
			destg[i]=app
		else
			destg[i]=db_tot[i,2]
	}
	for(i in 1:dim(db_tot)[[1]])
	{
		if(length(strsplit(db_tot[1,3],"->")[[1]])==2)
			srcdestg[i]<-paste(srcg[i],destg[i],sep="->")
		else if(length(strsplit(db_tot[1,3],"<->")[[1]])==2)
			srcdestg[i]<-paste(srcg[i],destg[i],sep="<->")
		else if(length(strsplit(db_tot[1,3],"-><-")[[1]])==2)
			srcdestg[i]<-paste(srcg[i],destg[i],sep="-><-")
	}
	out<-cbind(srcg,destg,srcdestg,labelg)
	out

}
inters_MI_string<-function(modid,string_info,mi_value)
{
	pos<-NULL
	modid1<-subset(modid, mi>mi_value)
	for(i in 1:length(modid1[,1]))
	{
		a1<-which(sub("G","",modid1[i,1])==string_info[,1])
		a2<-which(sub("G","",modid1[i,1])==string_info[,2])
		b1<-which(sub("G","",modid1[i,3])==string_info[,1])
		b2<-which(sub("G","",modid1[i,3])==string_info[,2])
		if(length(intersect(a1,b2))!=0 || length(intersect(a2,b2))!=0)
			pos<-c(pos,i)
	}
	return(modid1[pos,])
}

getedge <- function(outp_fc_gamma1,pathinfo) {
 edges00 <- src_dest(as(igraph.to.graphNEL(outp_fc_gamma1),"matrix"))
 edges01 <- sub("ENTREZID:","",edges00)
 edges02 <- sub("ENTREZID:","",edges01)
 edges0<-assignlabel(pathinfo,edges02) 
 return(edges0)
}


originalmodel <- function(edges0)
{
  mod0 <- ""
  for (i in 1:nrow(edges0)) {
    if (grepl("<->",edges0[1,3])) {
      lavmod <-  paste("G",edges0[i,1],"~~","G",edges0[i,2],sep = '')
    } else  {
      lavmod <-  paste("G",edges0[i,2],"~","G",edges0[i,1],sep = '')
    }
    
    mod0 <- paste(mod0,lavmod,sep=' \n ' )
  }
  return(mod0)
}

modifiedmodel <- function(mod0,MIchoice0)
{
  mod1 <- mod0
  for (i in 1:nrow(MIchoice0)) {
    lavmod <- paste(MIchoice0[i,1],MIchoice0[i,2],MIchoice0[i,3],sep = '')
    mod1 <- paste(mod1,lavmod,sep=" \n ")
  }
  return(mod1)
}

nodemodel <- function(modfinal,data_pathway)
{
  mod1g <- modfinal
  genenames <- colnames(data_pathway)[2:length(data_pathway)]
  for (i in 1:length(genenames)) {
    glavaan <- paste(genenames[i],'~group',sep = '')
    mod1g <- paste(mod1g,glavaan,sep=" \n ")
  }
  return(mod1g)
}

edgemodel <- function(edges0,MIchoice)
{
  mod0diff <- ""
  mod1diff <- ""
  for (i in 1:nrow(edges0)) {
    if (grepl("~~",edges0[i,3])) {
      lavmod0 <-  paste(edges0[i,1],"_0~~b",i,"*",edges0[i,2],"_0",sep = '')
      lavmod1 <-  paste(edges0[i,1],"_1~~a",i,"*",edges0[i,2],"_1",sep = '')
    } else  {
      lavmod0 <-  paste(edges0[i,1],"_0~b",i,"*",edges0[i,2],"_0",sep = '')
      lavmod1 <-  paste(edges0[i,1],"_1~a",i,"*",edges0[i,2],"_1", sep = '')
    }
    
    mod0diff <- paste(mod0diff,lavmod0,sep=' \n ' )
    mod1diff <- paste(mod1diff,lavmod1,sep=' \n ' )
  }
  
  if (!is.null(MIchoice)) {
    for (i in 1:length(MIchoice)) {
      MIpathway <- stringr::str_split(MIchoice[i], " ")[[1]]
      lavmod0 <- paste(MIpathway[1],"_0",MIpathway[2],"b",i+nrow(edges0),"*",MIpathway[3],"_0",sep = '')
      lavmod1 <- paste(MIpathway[1],"_1",MIpathway[2],"a",i+nrow(edges0),"*",MIpathway[3],"_1",sep = '')
      
      mod0diff <- paste(mod0diff,lavmod0,sep=' \n ' )
      mod1diff <- paste(mod1diff,lavmod1,sep=' \n ' )
    }
    
    moddiff <- paste(mod1diff,mod0diff,sep=' \n ' )
    for (i in 1:(nrow(edges0)+length(MIchoice))) {
      lavmod <- paste('dif',i,':=','a',i,'-b',i,sep = '')
      moddiff <- paste(moddiff,lavmod,sep=' \n ' )
    }
  }else{
    moddiff <- paste(mod1diff,mod0diff,sep=' \n ' )
    for (i in 1:nrow(edges0)) {
      lavmod <- paste('dif',i,':=','a',i,'-b',i,sep = '')
      moddiff <- paste(moddiff,lavmod,sep=' \n ' )
    }
  }
  
  return(moddiff)
}

ggmmodel <- function(edges0,MIchoice)
{
  
  modbg <- ""
  moddg <- ""
  for (i in 1:nrow(edges0)) {
    if (grepl("<->",edges0[1,3])) {
      modbgtemp <-  paste("G",edges0[i,1],"*","G",edges0[i,2],sep = '')
      modbg <- paste(modbg,modbgtemp,sep='+ \n ' )
    } else  {
      moddgtemp <-  paste("G",edges0[i,2],"~","G",edges0[i,1],sep = '')
      
      if (i==1) {
        moddg <- moddgtemp
        
      } else {
        moddg <- paste(moddg,moddgtemp,sep=', \n ' )
      }
    }
    
  }
  
  
  for (i in 1:nrow(MIchoice)) {
    modbgtemp <- paste(MIchoice[i,1],'~',MIchoice0[i,3],sep = '')
    moddg <- paste(moddg,modbgtemp,sep=", \n ")
  }
  if (modbg=='') {
    ggmcode <- paste('makeMG(dg=DG(',moddg,'))')
  } else {
    ggmcode <- paste('makeMG(dg=DG(',moddg,'), \n bg=UG( \n ~',modbg,'))')
  }
  
  return(ggmcode)
}

## Modify modi
modindices_new <- function(mod, ...) {
  mi_tbl <- lavaan::modindices(mod, ...)
  mi_tbl$lhs <- str_replace(mi_tbl$lhs, "^z", "")
  mi_tbl$rhs <- str_replace(mi_tbl$rhs, "^z", "")
  return(mi_tbl)
}
