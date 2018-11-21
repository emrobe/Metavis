from os.path import join, dirname
import datetime, scipy, json, re
import numpy as np
import pandas as pd
import cProfile

from bokeh.io import curdoc
from bokeh.layouts import row, column, widgetbox, Spacer
from bokeh.models import ColumnDataSource, DataRange1d, Select, CustomJS
from bokeh.palettes import Blues4, Set1
from bokeh.plotting import figure
from bokeh.models.widgets import NumberFormatter, CheckboxButtonGroup, Slider, Div, PreText, Tabs, Panel, DataTable, TableColumn, TextInput, Button, Paragraph
from bokeh.tile_providers import CARTODBPOSITRON

def update_hist(attr, old, new):
    active = new
    layout.children[1].children[0].tabs[0] = make_plot_hist(title, stats, active)

def update_kde(attr, old, new):
    active = new
    layout.children[1].children[0].tabs[1] = make_plot_kde(title, stats, active)

def make_plot_hist(title, stats, active=[0]):

    readlegend = "Reads ({:,})".format(stats['readstotal'])
    contiglegend = "Contigs ({:,})".format(stats['contigstotal'])
    genesfulllegend = "Genes, complete ({:,})".format(stats['genesfullstotal'])
    genespartiallegend = "Genes, truncated ({:,})".format(stats['genespartialtotal'])

    xaxisrange = 1000
    yaxistype = "log"
    yaxislabel = "something"

    p = figure(width=800, title=title, tools="pan,wheel_zoom,save,reset", y_axis_type="log", x_minor_ticks=10, x_range=(0,1000), y_range=(0.5, 10**7), background_fill_color="#f4f4f4", toolbar_location="below")


    # Test which plots should be plotted
    if 0 in active:
        # Draw Reads plot
        p.quad(top=stats['readhist'], bottom=0.001, left=stats['readedges'][:-1], right=stats['readedges'][1:],fill_color=Set1[0], line_color="white", alpha=0.5, legend=readlegend)
    if 1 in active:
        p.quad(top=stats['contighist'], bottom=0.001, left=stats['contigedges'][:-1], right=stats['contigedges'][1:],fill_color=Set1[1], line_color="white", alpha=0.5, legend=contiglegend)
    if 2 in active:
        p.quad(top=stats['genesfullhist'], bottom=0.001, left=stats['genesfulledges'][:-1], right=stats['genesfulledges'][1:],fill_color=Set1[5], line_color="white", alpha=0.5, legend=genesfulllegend)

    if 3 in active:
        p.quad(top=stats['genespartialhist'], bottom=0.001, left=stats['genespartialedges'][:-1], right=stats['genespartialedges'][1:],fill_color=Set1[4], line_color="white", alpha=0.5, legend=genespartiallegend)

    p.legend.location = "top_right"
    p.legend.background_fill_color = "#fefefe"
    p.xaxis.axis_label = 'Length (bp)'
    p.yaxis.axis_label = 'Count (log10)'
    p.grid.grid_line_color="white"
    p = Panel(child=p, title="Histogram")
    return p

def make_plot_kde(title, stats, active=[0]):

    readlegend = "Reads ({:,})".format(stats['readstotal'])
    contiglegend = "Contigs ({:,})".format(stats['contigstotal'])
    genesfulllegend = "Genes, complete ({:,})".format(stats['genesfullstotal'])
    #genespartiallegend = "Genes, truncated ({:,})".format(stats['genespartialtotal'])
    maxy = float(max(stats['readkdey']))
    p = figure(width=800, title=title, tools="pan,wheel_zoom,save,reset", x_minor_ticks=10, x_range=(0,10000), y_range=(0, maxy), background_fill_color="#f4f4f4", toolbar_location="below")

    # Test which plots should be plotted
    if 0 in active:
        # Draw Reads plot
        p.line(stats['readkdex'], stats['readkdey'], legend=readlegend)
        #p.quad(top=stats['readhist'], bottom=0.001, left=stats['readedges'][:-1], right=stats['readedges'][1:],fill_color=Set1[0], line_color="white", alpha=0.5, legend=readlegend)
    if 1 in active:
        p.line(stats['contigkdex'], stats['contigkdey'], legend=contiglegend)
        #p.quad(top=stats['contighist'], bottom=0.001, left=stats['contigedges'][:-1], right=stats['contigedges'][1:],fill_color=Set1[1], line_color="white", alpha=0.5, legend=contiglegend)
    if 2 in active:
        p.line(stats['genekdex'], stats['genekdey'], legend=genesfulllegend)
        #p.quad(top=stats['genesfullhist'], bottom=0.001, left=stats['genesfulledges'][:-1], right=stats['genesfulledges'][1:],fill_color=Set1[5], line_color="white", alpha=0.5, legend=genesfulllegend)

    if 3 in active:
        pass
        #p.quad(top=stats['genespartialhist'], bottom=0.001, left=stats['genespartialedges'][:-1], right=stats['genespartialedges'][1:],fill_color=Set1[4], line_color="white", alpha=0.5, legend=genespartiallegend)

    p.legend.location = "top_right"
    p.legend.background_fill_color = "#fefefe"
    p.xaxis.axis_label = 'Length (bp)'
    p.yaxis.axis_label = 'Count (log10)'
    p.grid.grid_line_color="white"
    p = Panel(child=p, title="Kernel density estimation")
    return p

def make_statistics_table(metrics):
    with open (join(dirname(__file__), 'data/metrics.j')) as filehandle:
        metrics = json.load(filehandle)
    source = ColumnDataSource(data=dict())
    source.data = metrics
    columns = [
        TableColumn(field="Metric", title="Metric"),
        TableColumn(field="Reads", title="Reads", formatter=NumberFormatter(format="0,0")),
        TableColumn(field="Contigs", title="Contigs", formatter=NumberFormatter(format="0,0")),
        TableColumn(field="Genes", title="Genes", formatter=NumberFormatter(format="0,0"))
    ]

    data_table = DataTable(source=source, columns=columns, width=600)
    return data_table

def update_annotation_exact(attr, old, new):
    # This function probablby needs a rewrite if too slow. if not new is an empty search box, else has input text
    if not new:
        # Display table header as PreText
        information_1 = "Showing {} of {} annotated genes total".format(maxindeces-1, totalrows)
        annotation_information_1 = PreText(text=information_1)
        funclayout.children[1] = annotation_information_1
        # Reset table to original
        annotation = ColumnDataSource(data=df.head(n=maxindeces))
        csv_download = ColumnDataSource(data=df)
        columns = [TableColumn(field=str(key), title=str(key)) for key in list(df)]
        annotation_table = DataTable(source=annotation, columns=columns, width=1800)
        funclayout.children[2] = annotation_table
        # Redo buttons to original
        funcbutton1 = Button(label="Download annotation (tsv)", button_type="primary")
        funcbutton1.callback = CustomJS(args=dict(source=csv_download.data),code=open(join(dirname(__file__), "download_tsv.js")).read())
        funcbutton2 = Button(label="Download sequences (Nuc, fasta)", button_type="success", disabled=True)
        funcbutton3 = Button(label="Download sequences (Prot, fasta)", button_type="success")
        funcbutton3.callback = CustomJS(args=dict(source=csv_download.data, seqs=seqs),code=open(join(dirname(__file__), "download_protseq.js")).read())
        funcbuttons = row(funcbutton1, funcbutton2, funcbutton3)
        funclayout.children[3] = funcbuttons

    else:
        # Display "Searching while searching"
        annotation_information_1 = PreText(text="Searching...(expect atleast 5 second search time pr. 10000 annotated genes)")
        funclayout.children[1] = annotation_information_1
        # Use template df pandas to search and create new pandas with search results
        current = df.where(df.apply(lambda row: row.astype(str).str.contains(new, case=False).any(), axis=1, result_type='broadcast')).dropna(how='all')
        hitcount= len(current)
        annotation = ColumnDataSource(data=current.head(n=maxindeces))
        csv_download = ColumnDataSource(data=current)
        columns = [TableColumn(field=str(key), title=str(key)) for key in list(current)]
        annotation_table = DataTable(source=annotation, columns=columns, width=1800)
        funclayout.children[2] = annotation_table
        # Display hit information
        tablecount = len(current.head(n=maxindeces))
        information_1 = "Showing {} of {} hits of {} annotated genes total".format(tablecount, hitcount,totalrows)
        annotation_information_1 = PreText(text=information_1)
        funclayout.children[1] = annotation_information_1
        # Change buttons
        funcbutton1 = Button(label="Download annotation of current hits (tsv)", button_type="primary")
        funcbutton1.callback = CustomJS(args=dict(source=csv_download.data),code=open(join(dirname(__file__), "download_tsv.js")).read())
        funcbutton2 = Button(label="Download sequences of current hits (Nuc, fasta)", button_type="success", disabled=True)
        funcbutton3 = Button(label="Download sequences of current hits (Prot, fasta)", button_type="success")
        funcbutton3.callback = CustomJS(args=dict(source=csv_download.data, seqs=seqs),code=open(join(dirname(__file__), "download_protseq.js")).read())
        funcbuttons = row(funcbutton1, funcbutton2, funcbutton3)
        funclayout.children[3] = funcbuttons

def update_annotation_regex(attr, old, new):
    pass

# Load data generated by the generate-script
with open (join(dirname(__file__), 'data/export.j')) as filehandle:
    stats = json.load(filehandle)


######### Sample Provenance and Meta data #########
info_text = Div(width=1800, height=800, text="""This tab will contain <ul><li>General sample information (name, owner, etc...)</li><li>Provenance data (Selected tools, versions, runtimes etc.)</li><li>Meta data (Date, geo, temp, etc.)</li><li>Other information</li></ul>""")

######### Sequence distribution #########
#with open (join(dirname(__file__), 'data/metrics.csv')) as filehandle:
#    metrics = json.load(filehandle)

# Set colors
Set1 = Set1[7]

# Create plots, set title
title = "Sequence distribution of sample: {}".format(stats['sampleid'])
p1 = make_plot_hist(title, stats)
p2 = make_plot_kde(title, stats)

# Spacer for correct alignments
leftmetrics = Spacer(width=500)
topwidgetbox = Spacer(width=400, height=50)
topmetrics = Spacer(width=600, height=70)

# First set of radiobuttons
text1 = PreText(text="View histogram for:")
buttons1 = CheckboxButtonGroup(name="Active plots:", labels=["Reads", "Contigs", "Genes, complete", "Genes, truncated"], active=[0])
buttons1.on_change('active',update_hist)

# Second set of radiobuttons
text2 = PreText(text="View KD-estimation for:")
buttons2 = CheckboxButtonGroup(name="Active plots:", labels=["Reads", "Contigs", "Genes, complete", "Genes, truncated"], active=[0])
buttons2.on_change('active',update_kde)

# Sliders and selection inputs
sliderx = Slider(title="X scale", value=0, start=-100, end=100, step=10)
yaxis = Select(title="Y axis type", value='log', options=['log','linear'])

# Gather controls in controls, put the two figures in tabs and create csv table
controls = column(topwidgetbox, widgetbox(text1, buttons1, yaxis, sliderx, text2, buttons2, width=400))
tabs = Tabs(tabs=[p1, p2])

# Metric table (right)
met = ""
metrictable = make_statistics_table(met)
metrics = column(topmetrics, metrictable, sizing_mode='fixed')

# Final layout grid
layout = row(controls, tabs, leftmetrics, metrics, sizing_mode='fixed')

### Taxonomic Classification
# Small hack, iframe needs height in absolute px in order to strecht vertically... ? 
tax_iframe = Div(width=1800, height=800, text="""<iframe src="https://s1.sfb.uit.no/public/data/krona.html" style="overflow:hidden;height:800px;width:100%" frameborder="0"> </iframe>""")
unfuckspacer1 = Spacer(height=800, width=1)
tax = row(unfuckspacer1, tax_iframe)
print(dirname(__file__))
### Functional Assignment
# Load data
with open (join(dirname(__file__), 'data/csvdescriptions.j')) as filehandle:
    df = json.load(filehandle)
# Create dataframe, transpose and inser rows (genes) as an additional column. This will serve as the main template
df = pd.DataFrame(df).T
# Make separate pandas for sequences, an remove from template to increase performance
seqs = {k: v for k,v in zip(df.index.tolist(), df['seq'].tolist())}
df = df.drop(['seq'], axis=1).dropna(how='all')
df.insert(0, 'Genes', df.index)
# Set max table entries and get number of columns
totalrows = len(df)
maxindeces = 1001
# Create different Bokeh DataSources from main template.
annotation = ColumnDataSource(data=df.head(n=maxindeces))
csv_download = ColumnDataSource(data=df)
# Create view table
columns = [TableColumn(field=str(key), title=str(key)) for key in list(df)]
annotation_table = DataTable(source=annotation, columns=columns, width=1800)
# Search box
annotation_search = TextInput(value="", title="Search for annotation: (case insensitive)")
annotation_search.on_change('value', update_annotation_exact)
# Table header text
information_1 = "Showing {} of {} annotated genes total".format(maxindeces-1, totalrows)
annotation_information_1 = PreText(text=information_1)
# Download buttons (Options: ‘default’, ‘primary’, ‘success’, ‘warning’, ‘danger’, ‘link’)
funcbutton1 = Button(label="Download annotation (tsv)", button_type="primary")
funcbutton1.callback = CustomJS(args=dict(source=csv_download.data),code=open(join(dirname(__file__), "download_tsv.js")).read())
funcbutton2 = Button(label="Download sequences (Nuc, fasta)", button_type="success", disabled=True)
funcbutton3 = Button(label="Download sequences (Prot, fasta)", button_type="success")
funcbutton3.callback = CustomJS(args=dict(source=csv_download.data, seqs=seqs),code=open(join(dirname(__file__), "download_protseq.js")).read())
funcbuttons = row(funcbutton1, funcbutton2, funcbutton3)
# Put it all together
funclayout = column(annotation_search, annotation_information_1, annotation_table, funcbuttons)


### Geographic Context
geo = figure(x_range=(-18000000, 18000000), y_range=(-7000000, 9000000),
           x_axis_type="mercator", y_axis_type="mercator")
geo.add_tile(CARTODBPOSITRON)

### Downloads 
files = ['sample_S001_r1.fastq.gz', 'sample_S002_r2.fastq.gz', 'sample_S001_r1_trimmed.fastq.gz', 'sample_S002_r2_trimmed.fastq.gz', '16S_predicted.fasta', 'final.contigs.fasta']
# Input reads
dlheader = row(Div(text="""<h2>In this section you can download all files relevant to your sample (temporarily disabled)</h2>""", width=1200))
dlcont1 = Div(text="<h3>Input Read R1: {}</h3>".format(files[0]))
dlbutton1 = Button(label="Download", button_type="success", disabled=True)
dlrow1 = row(dlcont1, dlbutton1)
dlcont2 = Div(text="<h3>Input Read R2: {}</h3>".format(files[1]))
dlbutton2 = Button(label="Download", button_type="success", disabled=True)
dlrow2 = row(dlcont2, dlbutton2)
# Trimmed reads
dlspacer1 = Spacer(height=50)
dlcont3 = Div(text="<h3>Trimmed  Read R1: {}</h3>".format(files[2]))
dlbutton3 = Button(label="Download", button_type="success", disabled=True)
dlrow3 = row(dlcont3, dlbutton3)
dlcont4 = Div(text="<h3>Trimmed  Read R2: {}</h3>".format(files[3]))
dlbutton4 = Button(label="Download", button_type="success", disabled=True)
dlrow4 = row(dlcont4, dlbutton4)
# rRNA
dlspacer2 = Spacer(height=50)
dlcont5 = Div(text="<h3>16S rRNA sequences: {}</h3>".format(files[4]))
dlbutton5 = Button(label="Download", button_type="success", disabled=True)
dlrow5 = row(dlcont5, dlbutton5)
# Contigs
dlspacer3 = Spacer(height=50)
dlcont6 = Div(text="<h3>Assembled contigs: {}</h3>".format(files[5]))
dlbutton6 = Button(label="Download", button_type="success", disabled=True)
dlrow6 = row(dlcont6, dlbutton6)
# Put all rows in one column
dl = column(dlheader, dlrow1, dlrow2, dlspacer1, dlrow3, dlrow4, dlspacer3, dlrow5, dlrow6)

### Make panels from all individial panel layouts
info_panel = Panel(child=info_text, title="Overview")
seq_panel = Panel(child=layout, title="Sequence Distribution")
tax_panel = Panel(child=tax, title="Taxonomic Classification")
func_panel = Panel(child=funclayout, title="Functional Assignment")
geo_panel = Panel(child=geo, title="Geographic Context")
dl_panel = Panel(child=dl, title="Downloads")

# Final layout grid
window = Tabs(tabs=[info_panel,seq_panel, tax_panel, func_panel, geo_panel, dl_panel])

# Hail mary
curdoc().add_root(window)
curdoc().title = "Sequence_distribution"
