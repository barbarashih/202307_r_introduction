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
# Download the data from
# https://raw.githubusercontent.com/barbarashih/202307_r_introduction/main/data/day2/gene_exprs_annotated.csv
plot_df <- read.delim("data/day3/gene_exprs_annotated.csv")
```
		
2. 
I want to make a plot where :
- x-axis indicates TERT expression 
- y-axis indicates age group 
- The size of the dot indicates age (bigger dot, older individual)
- The colour indicates SMTSD
It would look like the plot below.
![image](https://github.com/barbarashih/202307_r_introduction/assets/8283969/91bd3ede-2bcd-4946-9185-87517c1216ef)

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
ggplot(data=plot_df, aes(x=AGE, y=TRET, colour= SMTSD, size=AGE)) +
	geom_point( alpha = 0.5, position = position_jitter(w = 0.2))  +
	theme_bw()
	
```


## Tasks
### Task 1.
<details><summary>Hint 1 </summary> 
Look at food_count in Day 2 Spot the error 3, for inspiration on how you might set a starting variable to keep your count.
</details>
<details><summary>Hint 2 </summary> 
You can loop through the length of the input vector, and add on each element of the vector to a variable you set at the beginning for each loop.
</details>
<details><summary>Hint 3 </summary> 
Consider food_count in Day 2 Spot the error 3. In each loop, 1 was added to food_count. Can you add each element of a vector instead?
</details>
Create your own sum function.

### Task 2.
<details><summary>Hint 1 </summary> 
Similar to Task 1 on how to make your own sum function, except you would need to divide by the length to get mean.
</details>
Create your own mean function.

### Task 3 (in progress.. check back later)
Use the data from Day 3 Section 3, plot:
- These genes: TERT,ELOVL3,KRT79 as y-axis
- Age group as x-axis
- Make boxplot (try googling how to use ggplot to make boxplot)
- Use grid.draw(g) to make the plot side by side for the two SMTSD (i.e. right hand side plot is for lower leg, and left hand side plot is for sun protected)
- Save each plot with a different file name