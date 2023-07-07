# Day 2 homework
## Spot the errors
Copy/paste these code, and try to fix them to make working code

1. 
<details><summary>Hint 1 </summary> 
Use head(plot_df) to look at the beginning of plot_df.
</details>
<details><summary>Hint 2 </summary> 
What kind of separator is it? Use ?read.delim to find out what the default separator is.
</details>
<details><summary>Hint 3 </summary> 
You can use sep="," to change the default separator (see Task 5.4 for examples on how to import a .csv file)
</details>
```r
	# Read in dataframe
	# Download the data from https://raw.githubusercontent.com/barbarashih/202307_r_introduction/main/data/day2/gene_exprs_annotated.csv
	plot_df <- read.delim("data/day3/gene_exprs_annotated.csv")
```
2. 
I want to make a plot where :
- x-axis indicates TERT expression 
- y-axis indicates anatomical site 
- The size of the dot indicates age group (bigger dot, older individual)
It would look like the plot below.
<details><summary>Hint 1 (error about TRET) </summary> 
Double check the spelling of the gene name indicated in the code and the gene name in my example plot. Typos are a common source of error when plotting genes.
</details>
<details><summary>Hint 2 </summary> 
plot_df$AGE is currently a character vector. You need to create a column that is a numeric vector
</details>
<details><summary>Hint 3 </summary> 
You can use strsplit to split a value
</details>
<details><summary>Hint 4 </summary> 
Once you have done strsplit(plot_df$AGE, "-"), take the first value of each character vector in each element of the list (see Day3 Task 1.2).
</details>
<details><summary>Hint 5 </summary> 
Remember to convert the values from characters to numerics (see Day 1, Task 4.2).
</details>
<details><summary>Hint 6 </summary> 
Remember to change the plot value from AGE for size to this newly generated column. Make sure you haven't changed x axis away from AGE; if x-axis becomes numeric, the position=position_jitter(w=0.2) will have a minial effect. This bit of the code introudce random spread of points that is 0.2 in width. 
</details>
```r
	# There is no error after the first + sign
	ggplot(data=plot_df, aes(x=AGE, y=TRET, colour= SMTSD, size=AGE)) + geom_point( alpha = 0.5, position = position_jitter(w = 0.2))  + theme_bw()
	
```


## Tasks
### Task 1.
<details><summary>Hint 1 </summary> 
Look at food_count in Day 2 Spot the error 3, for inspiration on how you might tackle this.
</details>
Create your own sum function.

### Task 2.
<details><summary>Hint 1 </summary> 
Similar to Task 1 on how to make your own sum function, except you would need to divide by the length to get mean.
</details>
Create your own mean function.

### Task 3.
- Save each figure you have created in the loop.
- Make 3 plots, with the 3rd plot using gender on x-axis
- Arrange the 3 plots in 2x2 (i.e. plot a square)

