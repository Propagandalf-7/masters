# A python script to count the words of a given text string.


text = """
In the first chapter, the models of this dissertation are presented and a dimensionless variational form for each is derived. The model problems for this dissertation are then formulated.

In the second chapter, two results are investigated. This first is on the existence and uniqueness of solutions for vibration problems. An article is discussed that proves the existence and uniqueness of solutions for a general vibration problem using semi-group theory and assuming four assumptions. All the models in this dissertation is a special case of this general problem, so the theory can be explained and applied to an example. The weak variational form for a cantilever Timoshenko beam is derived. In this derivation, we also give the definitions of the spaces required. The theory is then applied to this example problem. We show that the four assumptions can be proven to be true for the cantilever Timoshenko beam model.

The next result is in Modal Analysis. Again the cantilever Timoshenko beam model is used to first explain the theory, before the general case is discussed. A trial solution is presented and substituted into the partial differential equations. This results into two problems, an eigenvalue problem and an ordinary differential equation. These two problems are solved and substitution can confirm then that the trial solutions are indeed solutions. The method obtains the same result as separation of variables, but avoids any division. Finally a formal series solution for the example is given. The general case is then discussed and follows the same idea. Modal analysis is crucial for this dissertation. For the comparisons later in this dissertation, a method is required to compare the different models. Modal analysis shows that is the eigenvalues and eigenfunctions of the models compare well, then the solutions of the models should compare well. And since the Fourier Coefficients gets small very fast, it is only required to look at the first few eigenvalues.

In the third chapter, two theoretical results for the Finite Element Method is discussed. The first result is on the convergence of the Galerkin Approximation. The general case contains various symbols rectangularly used in application, so an example of the cantilever Timoshenko beam model is used and the Galerkin Approximation is derived. Then the general case is discussed from an article. The results from the article are presented in a reduced format to avoid introducing too many new symbols. We also present the results in four new Theorems, to present the results of the article in a concise manner.

The second result is on the convergence of the eigenvalues and eigenfunctions when using the Finite Element Method. The results are presented in textbook. Here we use the notation of the dissertation to present the results. The results are also expanded in an attempt to better explain the theory.
"""

word_count = len(text.split())


print("The text contains {} words".format(word_count))