\# Sports Scheduling Optimization Using Integer Programming



This project addresses the problem of designing a double round-robin sports

schedule using integer linear programming techniques. The case study is based

on the South American FIFA World Cup 2018 qualifiers organized by CONMEBOL.



The objective is to generate fair and balanced fixtures by minimizing the

number of \*double breaks\*, defined as repeated home or away games within the

same FIFA window. An additional constraint is introduced to prevent teams from

playing against Argentina and Brazil in consecutive matches. Optimized fixtures

for both scenarios, with and without the Brazil–Argentina constraint, are

provided in the `figures/` directory.





\## Methodology

\- Integer linear programming formulation

\- Binary decision variables and linear constraints

\- Optimization using the SCIP solver

\- Evaluation of different scheduling schemes



\## Results

The proposed models are able to generate optimal schedules with zero double

breaks in short computational times. Several scheduling strategies (French,

English and Inverted schemes) achieve optimal solutions under the proposed

constraints.



\## Repository Structure

\- `data/`: input data describing teams and match constraints

\- `figures/`: visual representations of generated schedules

\- `analysis.ipynb`: Jupyter notebook containing the integer programming model and optimization experiments

\- `report/`: full academic report (in Spanish)



\## Inspiration



This project is inspired by the paper:



\*Scheduling the South American Qualifiers to the 2018 FIFA World Cup by integer programming\*  

Guillermo Durán, Mario Guajardo, Denis Sauré  



The original paper can be found at:  

https://www.sciencedirect.com/science/article/abs/pii/S0377221717303909



