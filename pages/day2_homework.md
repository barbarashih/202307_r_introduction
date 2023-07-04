# Day 2 homework
## Spot the errors
Copy/paste these code, and try to fix them to make working code

1. 
```r
    # Loop through a vector and print each value
    items <- c(4, 5, 1)
    for(i in items){
        current_item <- items[i]
        print(current_item)
    )

```
2. 
```r
    # Make a list
    myMeals = list(breakfast=c("toast", "egg" "coffee"), 
                    lunch=c("sandwich", "crisps"), 
                    dinner=c("pie", "apple"))
```
3. 
```r
    # Go through each element of the list and print out its content with the index number
    # i.e. The expected output looks like this
    # "breakfast"
    # "0  toast"
    # "1  egg"
    # "2  coffee"
    # "lunch"
    # "0  sandwich"
    # "1  crisps"

    # Loop through the meals
    for(idx in 1:length(myMeals)){
        current_meal <- names(myMeals)[idx]
        # Print meal 
        print(current_meal)
        # Go through each food within a meal and print out each food
        food_count <- 1
        for(food in myMeals[current_meal]){
            out_line <- paste0(food_count, "  " food)
            print(out_line)       
            food_count <- food_count + 1
        }
    }
```
4.
```r
    # Make a function that multiply two input values. The function will print out a message sayign what the answer is, as well as return the answer so you can save it as an object for further operations
    # Use the function to create 2 objects and add them together 
    # You should get 13 as the final value
    myFun_numbers <- function(x, y){
        out_val <- x*y
        print(paste0("Answer is ", out_val))
    }
    num1 <- myFun_numbers(2,2)
    num2 <- myFun_numbers(3,3)
    num1 + num2
    
```


## Tasks
### Task 1. Download data from GTEx
Go to https://gtexportal.org/home/datasets. 
1. Under GTEx Analysis V8 > RNA-Seq data, click on "Gene TPM by tissue"
1. Download gene_tpm_2017-06-05_v8_skin_sun_exposed_lower_leg.gct.gz
1. Extract the downloaded gtc.gz file (you can unzip .gz file using 7-zip)
1. Import the extracted file into R (you might want to either edit this .gct file, or ignore the first 2 lines when importing the file)
