# Reproducible research: version control and R

Questions 1-3 can be found in my logistic growth repo: https://github.com/jamesimcculloch/logistic_growth

### Question 4
"A script for simulating a random walk is provided in the question-4-code folder of this repo. Execute the code to produce the paths of two random walks. What do you observe?"

Each graph shows the path of a particle in space (according to its x and y coordinates) over time, with the time denoted by the blue shade of the line. There is limited similarity between the two random walks beyond the fact that they are both highly irregular and that they begin at the coordinates (0,0)! There is no preferential direction of movement of the particle, with it equally likely to move in any direction at any point. The point to which it moves is independent of the direction in which it had travelled to get to the point from which it is moving. The random walk therefore has the potential to take the form of a relative linear movement across space (with short random detours) or a cluster of movement whereby the particle by chance stays in the same region of space. 

"Investigate the term **random seeds**. What is a random seed and how does it work?"

A random seed is the starting point for a pseudo-random number generator. While the numbers which follow the original seed (the input value) are random, when the same seed is used again, you will see the same numbers in the output. This allows the simulation of a random process, such as Brownian motion, to be reproducible. In the code in the question-4-code file, the runif() function selects a totally random number within a certain range, which in this case is between 0 and 2π. This produces the angle at which the simulated particle moves at each time point. However, this means that the code is not reproducible, and each time it is executed it will produce a different output. Instead, we could change the code to require a random seed when choosing the (pseudo-)random number for the angle of motion. Then when the same random seed is input, the output random walk will be the same each time that code is executed. You can do this by including the function set.seed() before the "for" loop, including any arbitrarily-chosen number as the argument. In my case, I have chosen the number 22 as my random seed, for no particular reason. Now, as long as I specify set.seed(22) in my function, I will get the same random walk as the output. I have edited the random_walk.R file to make it reproducible under the commit message "random_walk.R made reproducible". 

The changes that I made to the "random_walk.R" to make it reproducible can be seen in the below image taken of the compare view of my commit with the parent code:

<img width="1440" alt="random_walk_compare" src="https://github.com/jamesimcculloch/reproducible-research_homework/assets/150149794/b4820a78-9c9b-4e04-84a3-00f118c562fa">

I have also uploaded this image file to the repo, with the file name "random_walk_compare.png"

### Question 5

"Import the data for dsDNA viruses taken from the Supplementary Materials of the original paper into Posit Cloud. How many rows and columns does the table have?"

I can use the following code to find the number of rows and the number of columns of the table:

```{r}
data <- read.csv("Cui_etal2014.csv")
nrow(data)
ncol(data)
```
The output was 33 rows and 13 columns. This is not including the header row. If one were to include the header row, one would have to initially specify header = F in the read.csv() function. Including the header, there are 34 rows. 

"What transformation can you use to fit a linear model to the data? Apply the transformation."

In order to fit a linear model to the data, the variables need to be normalised. In the paper, this was done by log-transforming both the virion volume and the genome length. I initially checked for normality in the untransformed data using the Shapiro-Wilk test (as n < 50), as below:

```{r}
shapiro.test(data$Virion.volume..nm.nm.nm.)
shapiro.test(data$Genome.length..kb.)
```
This gave p-values of 4.254e-09 and 1.213e-10 indicating that the distribution of both variables is statistically significantly different from normal. 

I then created new columns with the log-transformed values of the virion volume and genome length columns:

```{r}
data$log.Virion.volume..nm.nm.nm. <- log(data$Virion.volume..nm.nm.nm.)
data$log.Genome.length..kb. <- log(data$Genome.length..kb.)
```
I then tested these with the Shapiro-Wilk test:

```{r}
shapiro.test(data$log.Virion.volume..nm.nm.nm.)
shapiro.test(data$log.Genome.length..kb.)
```

This gave p-values of 0.0001813 and 0.2731 respectively. This means that the genome length data is now normal (or at least not significantly deviating from normal), but the virion volume data still deviates from normal. Other transformations were tested (reciprocal, squaring, cubing, square-rooting), but none provided a larger p-value. Therefore I kept the log transform and accepted that it might not be possible to better normalise the data given the small sample size. 

"Find the exponent (𝛼) and scaling factor (β) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in Table 2 of the paper, did you find the same values?"

I fitted the linear model as so:

```{r}
model <- lm(log.Virion.volume..nm.nm.nm. ~ log.Genome.length..kb., data = data)
summary(model)
```
This gave the following output:

```{r}
Call:
lm(formula = log.Virion.volume..nm.nm.nm. ~ log.Genome.length..kb., 
    data = data)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.8523 -1.2530 -0.1026  1.0739  2.0193 

Coefficients:
                       Estimate Std. Error t value Pr(>|t|)    
(Intercept)              7.0748     0.7693   9.196 2.28e-10 ***
log.Genome.length..kb.   1.5152     0.1725   8.784 6.44e-10 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.263 on 31 degrees of freedom
Multiple R-squared:  0.7134,	Adjusted R-squared:  0.7042 
F-statistic: 77.16 on 1 and 31 DF,  p-value: 6.438e-10
```
As can be seen from the coefficients, the slope is 1.5152 (1.52 to 3 significant figures), which denotes 𝛼. β is the intercept, which in the coefficients is given as 7.0748. This can be inverse-logged as so:

```{r}
exp(7.0748)
```
To give 1181.807, our value for β. The p-values for each are 2.28e-10 (β) and 6.44e-10 (𝛼), which are indeed statistically significant. Both 𝛼 and β are the same as given in the dsDNA row of Table 2 in the paper (denoted as the allometric exponent and scaling factor, respectively), which is promising! 

"Write the code to reproduce the figure shown."

I used the following code to reproduce the figure shown in the question booklet:

```{r}
install.packages("ggplot2")
library(ggplot2) #ensuring ggplot is loaded
ggplot(aes(x = log(Genome.length..kb.), y = log(Virion.volume..nm.nm.nm.)), data = data) +
  xlab("log [Genome length (kb)]") +
  ylab("log [Virion volume (nm3)]") + #x and y axis labels
  geom_point() + #adding points to the graph
  geom_smooth(method = "lm") + #adding the fitted linear model line to the graph
  theme_light() + #changing the theme to match that of the graph in the booklet
  theme(axis.title=element_text(face="bold")) #making the axis labels bold
```
I wasn't entirely sure on which theme was used in the booklet. I eventually settled on theme_light() as the one that looked most similar to me. But theme_bw() and theme_linedraw() also look similar. I would have applied these in the same way as I did with theme_light(). The size of my points is a little smaller, but this depends on how zoomed in the plot is; I took the screenshot from the 'zoom' view, which makes the points appear smaller relative to the rest of the figure. My reproduced version of the plot is below:

<img width="994" alt="reproduced_figure" src="https://github.com/jamesimcculloch/reproducible-research_homework/assets/150149794/3389f09a-3c4d-4e53-a382-37266de21ea5">

"What is the estimated volume of a 300 kb dsDNA virus?"

The following code was used to predict the virion volume of a 300 kb dsDNA virus:

```{r}
prediction <- predict(model, data.frame(log.Genome.length..kb. = log(300)))
exp(prediction)
```
This gave the following output:

```{r}
      1 
6698076 
```
Therefore the estimate volume of a 300 kb dsDNA virus is 6698076 nm<sup>3</sup>.

### Bonus question

"Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? What limitations do they have?"

Reproducibility and replicability are often confused or conflated. Reproducibility is whether the use of the same data and methods reaches the same results. Replicability is whether studies attempting to address the same question reach the same conclusions, but using different data. 

Git is a software system for version control, tracking changes in code. These changes can be collected in repositories, which can be made available on the website GitHub, alongside the code itself. This allows other scientists to access the code that you have used in a particular scientific investigation. Your data can also be uploaded into a GitHub repository and made publicly available. Git and GitHub can be used to enhance reproducibility because it allows other workers to use your data and code and attempt to produce the same results as in your analysis. 

In terms of replicability, workers in the same field attempting to answer the same question are able to use your code (as accessed via GitHub) to help with their analysis, as it could give inspiration as to how best to fit a certain model, carry out a particular statistical test, visualise the results, _et cetera_. Even the exact same code could be applied to different data, to properly test the validity of your methods and the replicability of your research. 

In addition to these benefits, the version control aspect of git means that methods can continuously be updated to reflect scientific advancements and the development of new computational techniques, packages, and methods generally, but the old versions of codes and methods are still available to test the reproducibility and replicability of investigations carried out prior to these advancements. 

Furthermore, GitHub provides the space for providing more detailed elaboration as to the methods used for a certain paper, as methods sections in papers are often required by the published to be brief. This is particularly true for the highest impact journals, the studies in which are often the most impactful and therefore most in need of being reproducible and replicable. 

The branching and merging features are especially good for collaboration. Branching allows other workers to use and test the same code for reproducibility, and modify it and apply it to different data to test replicability. Merging allows updated methods, reflecting scientific advancements, to be added to a base branch and become the new template for further investigation using the same/similar methods. 

The online repositories are backed up; if the code for a specific investigation were only stored on a single device, which was then compromised, it would be difficult or impossible to test the reproducibility or replicability of that investigation. Hosting of the repositories on the internet means that they can be accessed from any device.

However, there are also limitations to using git and GitHub. The actions on the website are not always intuitive for new users, and therefore uptake of GitHub by new scientists may not be easy. This is even worse when one is using git without the GitHub interface. In addition, the hosting of the repositories online make them vulnerable to security breaches, even when the settings of the repository mean that it is not generally publicly viewable. The online nature of GitHub also necessitates an internet connection, so it may not be possible to use this interface when on fieldwork in a remote location without internet connection. One could use git without the GitHub interface, but this requires more technical understanding.

## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   - A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
   - Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
   - Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
   - Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \beta L^{\alpha}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   - Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
   - What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
   - Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
   - Write the code to reproduce the figure shown below. (10 points)

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  - What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
