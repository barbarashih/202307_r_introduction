# Day 1 homework
## Spot the errors
Copy/paste these code. They each have one or more intentional mistakes. Try and fix them to make the code work

1. 
```r
    # Assign 12, 13 and 15 to three objects
    # add them together to make a forth object
    # and then print the final object
    a_num <- 12
    b num <- 13
    3c_num <- 15
    d.num <- a + b + 3c_num

    print(d_num)

```
2. 
```r
    # Download the file from here https://github.com/barbarashih/202307_r_introduction/blob/main/data/day1/homework_gene_expression.csv
    # Read it in
    gene_exprs_df <-- read.csv("C:\Users\shihb\OneDrive - Lancaster University\work\teaching\workshop\_github\r_introduction\data\homework_gene_expression.csv")

    # check what the dataframe looks like
    str(gene_exprs_df)
    head(gene_exprs_df]

    # Get the average for Brain3, Blood3 and Blood2
    # You can use ?mean to look at the options for the function
    mean(gene_exprs_df$Brain2)
    mean(gene_exprs_df$Blood3)
    mean(gene_exprs_df$Blood2)

    # Get the standard deviation for Brain2, Blood3 and Blood2
    # You can use ?sd to look at the options for the function
    sd(gene_exprs_df$Brain2)
    sd(gene_exprs_df$Blood3)
    sd(gene_exprs_df$Blood2)

```

## Tasks
#### Download and filter sample annotation data from GTEx project
Go to https://gtexportal.org/home/datasets. Under GTEx Analysis V8 > Annotations, download:
1) GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx (explainations for column names for GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)
2) GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt (the file is tab deliminted, so you would want to use sep="\t")

Import 2) into R Studio as a dataframe. 


Filter the dataframe so you only keep "Skin - Sun Exposed (Lower leg)" in the SMTSD column
<details><summary>Hint 1 </summary> 
Review the conditional operator section of the tutorial, where we used > to keep rows where Brain gene expression is higher than 0. 
</details>
<details><summary>Hint 2 </summary> 
In this case, use == for equal to compare the annotation columns. Look at the GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx to check which annotation column you would need to filter by.
</details>

Filter the dataframe so you keep only "RNASEQ" samples
<details><summary>Hint 3 </summary> 
Look at the GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx to check which annotation column you would need to filter by.
</details>


Save the filtered dataframe