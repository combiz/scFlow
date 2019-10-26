---
title: "<b>scflow</b> - Quality Control Report"
date: "26 October, 2019"
output:
  html_document:
    theme: "flatly"
    toc: false
    fig_caption: true
    keep_md: true
    includes:
      after_body: footer.html
params:
  metadata_path: false
---





## Sample metadata

A total of 20 metadata variables were imported from the sample sheet for this sample: -

<!--html_preserve--><div id="htmlwidget-792a41ff1b3150366150" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-792a41ff1b3150366150">{"x":{"filter":"none","data":[["individual","batch","group","diagnosis","coverage","batches","sex","age","PMI","duration","capdate","prepdate","seqdate","nucleicount","cdnaconc","libraryconc","MLS","RIN","ap","aplevel"],["<i>factor<\/i>","<i>factor<\/i>","<i>factor<\/i>","<i>factor<\/i>","<i>factor<\/i>","<i>integer<\/i>","<i>factor<\/i>","<i>integer<\/i>","<i>integer<\/i>","<i>integer<\/i>","<i>factor<\/i>","<i>factor<\/i>","<i>factor<\/i>","<i>numeric<\/i>","<i>numeric<\/i>","<i>numeric<\/i>","<i>integer<\/i>","<i>numeric<\/i>","<i>factor<\/i>","<i>factor<\/i>"],["MS542","2","High","MS","Low","1","F","76","12","35","20181030","20181031","201811","1.96","0.584","13.1","394","5.3","Posterior","2"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>metadata<\/th>\n      <th>class<\/th>\n      <th>values<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":5,"scrollX":false,"columnDefs":[{"className":"dt-left","targets":"_all"}],"order":[],"autoWidth":false,"orderClasses":false,"lengthMenu":[5,10,25,50,100]}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## Number of counts / features per cellular barcode
<div class = "row">
<div class = "col-md-6">
<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-2-1.png" alt="&lt;b&gt;Figure: Histogram of count depth per cell.&lt;/b&gt; A lower-limit threshold of 300 was applied (red line). "  />
<p class="caption"><b>Figure: Histogram of count depth per cell.</b> A lower-limit threshold of 300 was applied (red line). </p>
</div>
</div>

<div class = "col-md-6">
<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-3-1.png" alt="&lt;b&gt;Figure: Histogram of number of genes per cell.&lt;/b&gt; A lower-limit threshold of 100 was applied (red line). "  />
<p class="caption"><b>Figure: Histogram of number of genes per cell.</b> A lower-limit threshold of 100 was applied (red line). </p>
</div>
</div>
</div>

## Count depth distribution by barcode rank (high to low counts)

<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-4-1.png" alt="&lt;b&gt;Figure: Barcode count depth rank plot.&lt;/b&gt; The 'elbow' indicates where count depth decreases rapidly (relative increase in background counts), and can be used to inform the count depth threshold.  The applied lower-limit counts threshold is indicated at 300 counts (red line). "  />
<p class="caption"><b>Figure: Barcode count depth rank plot.</b> The 'elbow' indicates where count depth decreases rapidly (relative increase in background counts), and can be used to inform the count depth threshold.  The applied lower-limit counts threshold is indicated at 300 counts (red line). </p>
</div>

## Number of genes versus count depth 

<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-5-1.png" alt="&lt;b&gt;Figure: Number of genes versus count depth coloured by relative mitochondrial counts.&lt;/b&gt; The count-depth threshold of 300 counts and the number of genes threshold of 100 genes are indicated with vertical and horizontal red lines, respectively. Cells with high mitochondrial counts are typically in cells with relatively lower count depth. Cells with fractional mitochondrial counts higher than 0.1 (i.e. 10.00%) were filtered."  />
<p class="caption"><b>Figure: Number of genes versus count depth coloured by relative mitochondrial counts.</b> The count-depth threshold of 300 counts and the number of genes threshold of 100 genes are indicated with vertical and horizontal red lines, respectively. Cells with high mitochondrial counts are typically in cells with relatively lower count depth. Cells with fractional mitochondrial counts higher than 0.1 (i.e. 10.00%) were filtered.</p>
</div>


## Fraction of mitochondrial / ribosomal counts
<div class = "row">
<div class = "col-md-6">
<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-6-1.png" alt="&lt;b&gt;Figure: Histogram of mitochondrial fraction per cell.&lt;/b&gt; A upper-threshold of 0.1 (i.e. 10.00%) maximum mitochondrial faction was applied (red line)."  />
<p class="caption"><b>Figure: Histogram of mitochondrial fraction per cell.</b> A upper-threshold of 0.1 (i.e. 10.00%) maximum mitochondrial faction was applied (red line).</p>
</div>
</div>

<div class = "col-md-6">
<div class="figure" style="text-align: center">
<img src="/home/ckhozoie/Documents/scflow/qc_files/figure-html/unnamed-chunk-7-1.png" alt="&lt;b&gt;Figure: Histogram of ribosomal fraction per cell.&lt;/b&gt; A upper-threshold of 1 (i.e. 100.00%) maximum ribosomal faction was applied (red line)."  />
<p class="caption"><b>Figure: Histogram of ribosomal fraction per cell.</b> A upper-threshold of 1 (i.e. 100.00%) maximum ribosomal faction was applied (red line).</p>
</div>
</div>
</div>
