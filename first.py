
import numpy as NP
import pickle
import matplotlib as MP
import pylab as PL




# First method to determine non scaling genes

key_array=NP.loadtxt('/Users/ilhamharmach/Desktop/FIRST_Python_projevt/Non_scaling.txt',dtype='str')

GFP=pickle.load(open('/Users/ilhamharmach/Desktop/FIRST_Python_projevt/mean_GFP_of_genes_BY.p','rb'))




values=[]

# values will contain the name of the gene and its corresponding expression
for i in key_array:
    values.append((i,GFP[i][0]))

# Sorting the genes depending on which have the highest expression
dtype =[('name','S10'),('gene_expression',float)]
a=NP.array(values,dtype=dtype)
b = NP.sort(a,order='gene_expression')
# c is our sorted list of genes who have the highest expression
c=sorted(b,reverse=True)


# Save my list in a txt file in order to have an other list


first_method_list= open('first_method.txt', 'w')

for item in c:
    first_method_list.write("%s\n" % item)


#  If I didn't want to have each gene in a different line
#  with open("file_1.txt", "w") as output:
#  output.write(str(c))




# Second METHOD :

# Part 1: I have a dictionary and I want to get the Y intercepts in a list

# My file is called Final fit dictionary , it contains the intercepts
# I want to create a list of all the intercepts

Final_fit=pickle.load(open('/Users/ilhamharmach/Desktop/FIRST_Python_projevt/final_fit_dic_BY.p','rb'))


# This is my list of intercepts
my_file=NP.loadtxt('/Users/ilhamharmach/Desktop/FIRST_Python_projevt/Yangetal_supp1.txt',skiprows=1,dtype='str')

# Turning my text file into an array to be able to use the "where" function

my_file=NP.array(my_file)

# Creating a list that will contain the name of the gene,its corresponding intercept and its 130

new_table=[]
list1=[]


# In order to plot the intercepts and the evolution of the cells after 170
# We search in all the genes we want, see if they coincide with the cells in my_file
for i in Final_fit.keys():
    a=NP.where(my_file==i)[0]
    if NP.size(a)!=0: # if the gene is actually present in my_file
        b=NP.where(my_file==i)[0][0]# I have the exact line of the gene that I need
        print b
        if my_file[b,2]!='NA':      #Get rid of all the missing data
            if abs(float(my_file[b,2])) < 4:
                list1.append(Final_fit[i][0][0])    # add the genes that correspond to all the previous conditions
                new_table.append(float(my_file[b,2])) # I have my list of 130
            


#plotting the intercepts and the 130

PL.plot(list1,new_table,'.','r')
)
PL.axis([-0.5,2,-3,3])
PL.xlabel("Non scaling metric Y intercepts")
PL.ylabel("Partitioning Z score /130 minutes")
PL.show()


# Programme pour plotter un histogramme avec les différentes valeurs des intercepts de la list1 des aanciens data;
# On se débarasse des baleurs qu'on aime pas et apres on le plotte, et on peut meme ploter une valeur d'intercept`#
# de nos nouveaux genes qu'on a avec les resultats du Flow cytometer
a=NP.array(list1)
k=a[NP.where(a<3)]
h=k[NP.where(k>-3)]
weights = NP.ones_like(h)/len(h)
PL.hist(h,100,weights=weights)
PL.show
# Pour plotter un intercept d'un gène mesuré sur l'histogramme
PL.plot(.86*[1,1],[0,.12],'r')
# Pour avoir le poucentage de gènes dont l'intercept est plus grand que 0.5
float(len(a[NP.where(a>.5)]))/float(len(a))



#Rajouter les axes et et les légéndes
# x=== non scaling metric Y intercepts
# y=Partitioning Z score /130 minutes

# In order to plot the intercepts and the evolution of the cells after 170
# We search in all the genes we want, see if they coincide with the cells in my_file
for i in Final_fit.keys():
    a=NP.where(my_file==i)[0]
    if NP.size(a)!=0: # if the gene is actually present in my_file
        b=NP.where(my_file==i)[0][0]# I have the exact line of the gene that I need
        print b
        if my_file[b,3]!='NA':      #Get rid of all the missing data
            list1.append(Final_fit[i][0][0])    # add the genes that correspond to all the previous conditions
            new_table.append(float(my_file[b,3])) # I have my list of 170
            
            


#Plotting the intercepts and the 170 results

PL.plot(list1,new_table,'.','r')
PL.axis([-0.5,2,-3,3])
PL.xlabel("Non scaling metric Y intercepts")
PL.ylabel("Partitioning Z score /170 minutes")
PL.show()




new_table=[]
list1=[]



# In order to plot the intercepts and the evolution of the cells after 170
# We search in all the genes we want, see if they coincide with the cells in my_file
for i in Final_fit.keys():
    a=NP.where(my_file==i)[0]
    if NP.size(a)!=0: # if the gene is actually present in my_file
        b=NP.where(my_file==i)[0][0]# I have the exact line of the gene that I need
        print b
        if my_file[b,4]!='NA':      #Get rid of all the missing data
            list1.append()    # add the genes that correspond to all the previous conditions
            new_table.append(float(my_file[b,4])) # I have my list of overnight



#plotting the intercepts and the overnight

PL.plot(list1,new_table,'.','r')
PL.axis([-0.5,2,-3,3])
PL.xlabel("Non scaling metric Y intercepts")
PL.ylabel("Partitioning Z score /overnight")
PL.show()

# if my_file[i,2]!='N/A': # It doesn't work because I use a name of a gene (i) and it is not a dictionary



# In order to plot the intercepts and the evolution of the cells after 170
# We search in all the genes we want, see if they coincide with the cells in my_file
for i in Final_fit.keys():
    a=NP.where(my_file==i)[0]
    if NP.size(a)!=0: # if the gene is actually present in my_file
        b=NP.where(my_file==i)[0][0]# I have the exact line of the gene that I need
        print b
        if my_file[b,3]!='NA':      #Get rid of all the missing data
            list1.append(Final_fit[i][0][0])    # add the genes that correspond to all the previous conditions
            new_table.append(float(my_file[b,3])) # I have my list of 170

                





list1=[]
new_table=[]
gene_expression_list=[]
d_list=[]
# In order to plot the intercepts and the evolution of the cells after 170
# We search in all the genes we want, see if they coincide with the cells in my_file
for i in key_array:
    a=NP.where(my_file==i)[0]
    if NP.size(a)!=0: # if the gene is actually present in my_file
        b=NP.where(my_file==i)[0][0]# I have the exact line of the gene that I need
        print b
        if my_file[b,3]!='NA':      #Get rid of all the missing data
            list1.append(i)    # add the genes that correspond to all the previous conditions
            new_table.append(float(my_file[b,3]))# I have my list of 170  for NON scaling genes
            gene_expression_list.append((i,GFP[i][0]))
            d_list.append((i,float(my_file[b,3]),GFP[i][0]))


# I have my whole list sorted and I can pick up my genes from the bottom
d.sort(key = operator.itemgetter(1, 2))




dtype =[('name','S10'),('Z-score',float),('gene_expression',float)]
e=NP.array(d_list,dtype=dtype)
f = NP.sort(e,order='Z-score')
g=sorted(f,reverse=True)






#gene_expression_list_sorted=sorted(gene_expression_list,reverse=True)



#new_list contains the name of the genes and their Zscore

new_list=[]
for i in list1:
    new_list.append((i,new_table[list1.index(i)]))

dtype =[('name','S10'),('Z-Score',float)]
a=NP.array(new_list,dtype=dtype)
b = NP.sort(a,order='Z-Score')

# c is our sorted list of genes who have the highest Zscore
c=sorted(b,reverse=True)
third_method_list= open('third_method.txt', 'w')

for item in c:
    third_method_list.write("%s\n" % item)

def column(matrix, i):
    return [row[i] for row in matrix]

a=[i for i in range(len(new_table))]

PL.plot(column(g,1),column(c,1),'.')
PL.axis([0,25000,-5,10])
PL.xlabel("Gene_expression")
PL.ylabel("Partitioning Z score for NON scaling genes /170")
PL.show()








def plot_fitted_graph(gene_name,lo_bound_hap,hi_bound_hap,file_type = 'pdf' ):
     final_data = pickle.load(open("all_data.p","rb"))
     sample=gene_name
     x=final_data[sample]['proc_data'][:,1]      # SSC
     y=final_data[sample]['proc_data'][:,3]      # Corrected GFP
     a=final_data[sample]['binned data'][:,0]    # binned data
     b=final_data[sample]['binned data'][:,3]
     PL.plot(x, y, '.', markersize = 12, alpha = .05, label = "{0}".format(sample))
     PL.plot(a,b'-k', linewidth = 2)
     ind_lo = NP.where(x> lo_bound_hap)[0]
     ind_hi = NP.where(x< hi_bound_hap)[0]
     cut_ind = NP.intersect1d(ind_lo,ind_hi)
     new_data_SSC_=x[cut_ind]
     new_data_GFP=y[cut_ind]

     print new_data_SSC.shape
     print new_data_GFP.shape


     slope, intercept =linear_model(new_data_SSC,new_data_GFP)
     PL.plot(new_data_SSC, intercept+slope*new_data_SSC, 'r')
     PL.legend(loc = "upper left")
     PL.xlabel("SSC (Volume)")
     PL.ylabel("GFP Intensity (AU)")
     PL.ylim(0,1050)
     PL.xlim(0,1050)
     savepath = '{0}.{1}'.format(sample,file_type)
     PL.savefig(savepath)
     PL.close()





# Plotting gene expression and Z score in the same plot and looking for genes that have the
# highest of both of them because they're most likely to be non scaling
# Starting with 2 of them
# We're taking :
# 'HMO1'
# 'YMR295C'
# 'PRY3'
# 'HHO1'

# "HM01" :YDR174W,H8
# YMR295C ,F7
#'PRY3',YJL078C, F6
# HHO1,  YPL127C, F11
# Going to the fridge and picking up frozen cells for which HMO1 is marked

# plate number
# SGD CORRESPONDING NAME
# well : H8 : compter les lettres verticalement et les chiffres horizontalement
# Faire ça pour mes 4 gènes













