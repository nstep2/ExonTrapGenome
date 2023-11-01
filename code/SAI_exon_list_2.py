



import scipy
import scipy.stats
res = scipy.stats.linregress(np.log10(x), np.array(y))

import statsmodels.formula.api as smf
import statsmodels.api as sm




outdir = exp_output_path.Compute_spliceAI_exon_id_scratch 
pdf_plots = PdfPages(outdir+'Compute_spliceAI_exon_id_scratch%d_reads.pdf' % (exon_count_build))









if "load" == "no":
    from experiment_paths.experiment_paths import *
    import pickle
    
    pickle_path = exp_output_path.SpliceAI_ET_exons + "SpliceAI_ET_exons.pickle"
    
    with open(pickle_path, "rb") as input_file:
        exon_id_list_spliceAI_scores_dict = pickle.load( input_file)
        










pickle_path = exp_output_path.SpliceAI_ET_exons + "SpliceAI_ET_exons.pickle"

with open(pickle_path, "rb") as input_file:
    exon_id_list_spliceAI_scores_dict = pickle.load(input_file)
    





for exon_id_list_name in exon_id_list_spliceAI_scores_dict:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[exon_id_list_name]  
    
    
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)



pdf_plots.close()








def get_bin_median(x,y,xlabel,ylabel,title,x_lower_lim,y_lower_lim,num_bins,pdf_plots):
    x_c = np.log10(x)
    y_c = np.log10(np.array(y) ) #*np.array(z))
    nx, ny = (num_bins, num_bins) #31
    xb = np.linspace(2, 5, nx)
    yb = np.linspace(y_lower_lim, 0, ny)
    
    
    
    
    bins_step = 0.25
    x_bins = np.arange(2,4.75,bins_step)
    y_vals = list()
    for bb in x_bins:
        y_vals.append(list())
        
    for ii, val in enumerate(x_c):
        for jj, x_val in enumerate(x_bins):
        
            if  x_c[ii] >= x_val and x_c[ii] < x_val + bins_step:
                y_vals[jj].append(y_c[ii])
    
    



    
    
    import scipy
    import scipy.stats
    res = scipy.stats.linregress(np.log10(x), np.array(y))
    
    
    import statsmodels.formula.api as smf
    import statsmodels.api as sm

    X = sm.add_constant(np.log10(x), prepend=False)
    Y=np.array(y)

    df = pd.DataFrame({'x':np.log10(x), 'y':np.array(y)})
    model =smf.ols(formula='y ~ x',data = df).fit()
    
    
    data_list = list()
    fig = plt.figure()
    for ii, y_list in enumerate(y_vals):
        median = np.median(y_list)
        quantiles = np.quantile(y_list,[0.25,0.75])
        spread = np.quantile(y_list,[0.1,0.9])
        bin_count = len(y_list)
        data_list.append([x_bins[ii],median,quantiles[0],quantiles[1],spread[0],spread[1], bin_count])
      
    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2, bin_count = zip(*data_list)

    
    return data_x,data_y, [res, model], bin_count







outdir = exp_output_path.Compute_spliceAI_exon_id_scratch 
pdf_plots = PdfPages(outdir+'combined_Compute_spliceAI_exon_id_scratch%d_reads.pdf' % (exon_count_build))


x,y,z = zip(*exon_id_list_spliceAI_scores_dict['Intronic'])




median_lines=dict()
median_lines['Intergenic']= dict()
median_lines['Intronic'] = dict()
median_lines['mRNA']     = dict()
median_lines['lncRNA']   = dict()
median_lines['Antisense']   = dict()


for exon_id_list_name in exon_id_list_spliceAI_scores_dict:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[exon_id_list_name]  
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)
            
    
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 3ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=25
    data_x,data_y,r_sq, bin_count = get_bin_median(x, np.array(y), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    median_lines[exon_id_list_name]['3ss'] = data_x,data_y
    median_lines[exon_id_list_name]['3ss_rsq'] = r_sq
    median_lines[exon_id_list_name]['3ss_bin_count'] = bin_count
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 5ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=25
    data_x,data_y,r_sq, bin_count = get_bin_median(x, np.array(z), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    median_lines[exon_id_list_name]['5ss'] = data_x,data_y
    median_lines[exon_id_list_name]['5ss_rsq'] = r_sq
    median_lines[exon_id_list_name]['5ss_bin_count'] = bin_count





print('\nStandard Error of the Regression for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.scale)) )


print('\nStandard Error of the Regression for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.scale)) )



print('\nrsq for OLS model for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.rsquared)) )

print('\nrsq for OLS model for 5ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.rsquared)) )


print('\nrsq for scipy linear fit for 3ss')
for exon_id_list_name in median_lines:
    res, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, res.rvalue**2) )


print('\nrsq for OLS model for 5ss')
for exon_id_list_name in median_lines:
    res, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, res.rvalue**2) )









fig,ax = plt.subplots()
for key in ['mRNA','lncRNA','Intergenic']:
    res, model = median_lines[key]['3ss_rsq']
    r_squared = res.rvalue**2
    data_x,data_y = median_lines[key]['3ss']
    plt.plot(data_x,10**np.array(data_y),label=key+" {:.2}".format(r_squared), color=color_dict[key])

         
plt.title('median 3ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 3\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)






fig,ax = plt.subplots()
for key in ['mRNA','lncRNA','Intergenic']:
    res, model = median_lines[key]['5ss_rsq']
    r_squared = res.rvalue**2
    data_x,data_y = median_lines[key]['5ss']
    plt.plot(data_x,10**np.array(data_y),label=key+" {:.2}".format(r_squared), color=color_dict[key])
#plt.plot([2,5],[(0.2),(0.2)], linestyle='dashed', alpha =0.5, color='orange',label='spliceAI 0.2 threshold')
plt.title('median 5ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 5\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)



pdf_plots.close()



data_x,data_y = median_lines['mRNA']['5ss']
print(  "Bins\t" + '\t'.join("{:,}".format(int(10**item)) for item in data_x) +'\n')

for key in ['mRNA','lncRNA','Intergenic']:
    bin_5ss = median_lines[key]['5ss_bin_count']
    print(key + '\t' + '\t'.join("{:,}".format((item)) for item in bin_5ss) +'\n')

    







for key in ['mRNA','lncRNA','Intergenic']:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[key]
    
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)
    
    
    
    with open( exp_output_path.Compute_spliceAI_exon_id_scratch  + "4G_{:}.txt".format(key), 'wt') as f:
        f.write('exon_count\tscores\n')

        for ii, val in enumerate(x):
            out_string = "{:}\t{:}\n".format(x[ii],y[ii])
            f.write(out_string)
    



for key in ['mRNA','lncRNA','Intergenic']:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[key]
    
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)
    
    
    
    with open( exp_output_path.Compute_spliceAI_exon_id_scratch  + "4H_{:}.txt".format(key), 'wt') as f:
        f.write('exon_count\tscores\n')

        for ii, val in enumerate(x):
            out_string = "{:}\t{:}\n".format(x[ii],z[ii])
            f.write(out_string)
    






    



