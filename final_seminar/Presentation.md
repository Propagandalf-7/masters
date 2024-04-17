

#########################################################################################################

Thank you everyone for joining my master's seminar today. This presentation is titled 
`Finite Element Analysis of multi-dimensional and simplified models for beams`. My thesis also covers 
plate models; however, today's presentation will focus solely on beams, despite the mention of plates 
in the seminar invitation, which was an error.

Beams are fundamental structural elements in engineering. They are in a vast array of applications 
including bridges, buildings, and even the International Space Station.

Simplified models are often used over more realistic models for vibration problems due to their 
reduced complexity and lower computational demands compared. More realistic models are also often 
higher dimensional models. However, the accuracy of these simplified models in practical applications
is not always certain. Today, we will examine some of the key factors that affect how well these models 
perform.

In this presentation, we will assess the validity of a Timoshenko cantilever beam model and a 
two-dimensional cantilever beam model. Here, 'validity' refers to how well the solution from the simplified 
model compares to the solution from a more complex, higher dimensional model. To be able to compare the 
models, we'll use modal analysis that will allow us to compare the eigenvalues and eigenfunctions. Modal 
analysis also tells us that if these eigenvalues compare well, then the solutions should also compare well.

The structure of this presentation is as follows:

    1. We will begin with a discussion on two literature studies. The first article, titled 
    `Comparison of Linear Beam Theories`, compares a cantilever Timoshenko beam with a two-dimensional beam.
    We will review the theory used by the authors and some of their results, and we will extend their work
    to a comparison of a two and three-dimensional model..

    2. Next, we will look at a study where researchers vibrated an actual beam and measured its natural 
    frequencies, comparing these to results to the Timoshenko beam theory and a three-dimensional beam 
    model using the Finite Element Methods.

    3. We will then quickly look at the models used in this presentation: the Timoshenko cantilever beam,
    the two-dimensional cantilever beam, and the three-dimensional cantilever beam.

    4. Following that, we will touch on some theory. We'll briefly discuss the well-posedness of the models
    without going into much detail. We will also then look at modal analysis, which is essential for comparing
    our models.
    
    5. Finally, we will present the numerical results. 
        - We will start by discussing the computation of eigenvalues and eigenfunctions using the 
        Finite Element Method.
        - We will then review the validity of the Timoshenko beam based on the findings from the
        first literature study, recalculating and verifying these results.
        - Finally we will extend our discussion to the validity of the two-dimensional beam model,
        following a similar methodology than the article.

Finally we will conclude our findings, and discuss future work.

#########################################################################################################

Conside the article `Comparison of linear beam theories` by Labuschagne, Van Rensburg and Van Der Merwe.

In this article, the authors investigate the efficacy of the Timoshenko beam theory against a more
realistic two-dimensional model. Specifically, the authors compare a cantilever Timoshenko beam to 
a cantilever two-dimensional beam, as a cantilever is use often in practical applications. 

The authors use modal analysis and the Finite Element Method to calculate the eigenvalues and
eigenfunctions of the models and compare these. We'll go into more detail into how the eigenvalues are
compared later in this presentation.

The authors show that the Timoshenko beam compares very well to the two-dimensional model. Especially for
the first ten eigenvalues. The authors also show that the shape of the models significantly influences 
the comparison. We'll also see this when we get to the numerical results.

The authors do not include damping in their investigation, but show that for practical applications, the
Timshenko beam model compares well to a two-dimensional model. Practical applications do refer to common
beam shapes. 

Although a direct comparison to a three-dimensional model is preferred, the authors recommend that the
two-dimensional model be used as an intermediate step, as this should mitigate some complexities that we'll
discuss in the numerical results.

In this article we will dicuss and replicate the results of the article, up to 5 significant digits. We
also extend the results to a comparison of a two-dimensional cantilever beam to a three-dimensional
cantilever beam.

#########################################################################################################

Next we look at an article we'll just quickly discuss that I included for interesting results. The article
is titled `On the valid frequency range of Timoshenko beam theory`. In this article, the authors report 
on an experiment that they conducted. They suspended a small beam, made from an aluminum alloy, at both 
endpoints on carbon fibre strings.

The beam is vibrated and a free-free beam is simulated. They then sweep through the frequencies, and 
measured the natural frequencies of the beam.

The results obtained is then compared to some theoretical results. 

Here is just a representation of the suspended beam.

And here we have the results from the article. The table shows the results of the first 8 eigenvalues. 
For each eigenvalue, they measured a `stiff` plane and a `flexible` plane, since the beam is not uniformly
square. ***SHOW BACK THE WIDTH AND BREADTH***

This table also then includes the results of the three-dimensional model using FEM, and the Timoshenko 
beam theory. And as we can see in the table, the results are very close. The three-dimensional model
does compare better, but the Timoshenko beam is also respectable.

#########################################################################################################

Now we can look at the models that we'll be using in this presentation. Note that we do not consider any
damping in the models, and also the models are all in a dimensionless form. We also will only look at 
a cantilever configuration with a sqaure cross-secion of the models, following the first article 
we discussed.

#########################################################################################################

We start with the Timoshenko cantilever beam.

Here we have the `equations of motion` and the `constitutive equations` for the Timoshenko beam model. 
In these equations, w is the displacement, phi is the angle of rotation, M is the bending moment, V is 
the shear force, Q is an external distributed load and alpha and beta are dimensionless constants.

Here is just a representation of the cantilever configuration and the boundary conditions. The
displacement w and rotation of the cross-secition phi is 0 at x =0, and the moment and shear force are 0
at x = 1.

#########################################################################################################

Next we have the two-dimensional cantilever model.

Here we have the `equations of motion` and the `constitutive equations`. Here u is the displacement vector.
T is the stress-tensor, Q is a distributed load, Epsilon is youngs modulus, nu is poissons ratio and gamma
is a dimensionless constant.

Here we also just have a representation of the two-dimensional cantilever beam with the boundary conditions.
In the boundary conditions, u is 0 at x_1 = 0 and the stress Tensor times the outward normal vector is 
0 on the rest of the body.

#########################################################################################################

Finally we have the three-dimensional cantilever model.

The `equations of motion` and the `constitutive equations` as well as the variables and constants are
similar to the two-dimensional model, so I will not go over them again, but I'll show them here on the 
screen quickly. 

The two-dimensional model is actually a special case of plane stress of the three-dimensional model. 
For this special case, we assume that all the stress tensor components in the 3rd axis is 0, the derivatives
in the direction of x3 is zero, as well as the strain component epsilon 33 is 0.

Is this is assumed, the two-dimensional model is obtained. This is not a general case of plane stress, as the
component epsion33 is not necessarily zero for the general case.

And then again we have a representation of the cantilever three-dimensional model. Again the boundary
conditions are similar to the two-dimensional case.

#########################################################################################################

Next we look at some theory for the models. For this theory, we'll consider a general vibration problem.

First we define the spaces where the general vibration problem is defined. This is from an article titled
Analysis of Solvability of Linear Vibration Models.

Let V, W, and X be real Hilbert spaces such that W is a linear subspace of X, and V is a linear subspace 
of W. X is called the global space, W is the inertia space and V is called the energy space. Each of these
spaces have an inner product and a induced norm.

Before the general vibration problem can be given, we need to define some bilinear forms. The bilinear forms
are a, b and c where a is the inner product of X, b is the inner product of V and c is the inner product of
W.

Then finally we can define our general vibration problem. Find a funcyion u in C such that u prime is 
continuous at 0 with respect to the norm of W, and for each value of t in an interval J, u(t) is in V, 
u'(t) is in V and u''(t) is in W, satisfying equation 2 and some initial conditions.

All of the models in this presentation are second order hyperbolic type problems, and each is a special
case of this general vibration problem. And since we assume no damping, we actually have that the 
bilinear form a is always 0. So we can write the general vibration problem without the term for the 
bilinear form a.

#########################################################################################################

Next we look at the well-posedness of the models. We wont go into too much detail, but rather just look 
at a result. In the article Analysis of the Solvability of Linear Vibration Models, the authors present 
results for the existence and uniqueness of solutions to the general vibration problem. 

The following assumptions are made in the article:

    A1 - V is dense in W and W is dense in X
    A2 - The norm of X is bounded above by the norm of W
    A3 - The norm of W is bounded above by the norm of V
    A4 - The bilinear form is non-negative, symmetric and bounded on V

Now since we have no damping and the bilinear form a is always 0, the forth assumption is automatically 
satisfied by our models.

The authors present the following theorem for the extence and uniqueness of solutions for the general 
vibration problem. For our models, it is therefore only required to prove that the assumptions are met 
and that the initial values are admissable. We will not look into this in this presentation.

#########################################################################################################

Next we will look into Modal Analysis. 

Consider the general vibration problem again. This time with the bilinear form a equal to 0 and no 
external force. 

Consider a trial solution for this problem. Let us take
u(t) = T(t)x with x in V and x not 0.

If we substitute this trial solution into our equation 4, we get the following. And due to the linearity
of the biliear form b, we can take out the T(t) on the left hand side and the T double prime (t) on the 
right hand side.

Dividing both sides by T(t) gives the following equation.

Form this it follows that this relation T''(t) over T(t) must be a constant.

Set this value equal to some number, let call it lambda. And at this point, the existence of such a value
lambda is not known. 

So we end up with two problems to solve. 
The first is called an eigenvalue problem. Here we need to find a real number lambda and a non zero x in
V such that equation 5 is satisfied.

We also have an ordinary differential equation for T(t) that we need to solve in equation 6.

#########################################################################################################

The article `Using energy methods to compare linear vibration models` by Civin, Van Rensburg and 
Van Der Merwe provides the solution to these two problems. In this article, the authors add an 
additional assumption, A5 to our previous assumptions. A5 states that the embedding of V into W is compact.

Using these assumptions, the authors prove that there exists a complete orthonormal sequence of eigenvectors
for the eigenvalue problem with a corresponding sequence of eigenvalues. These eigenvalues are positive and 
the orthogonality is with respect to the bilinear form c. Also the sequence of normalized eigenvectors 
x_i forms an orthonormal basis in W and sequence of eigenvalues {\lambda_i} is an infinite sequence 
with lambda_n tends to infinity as n tends to infinity. 
            
Following these results, for any u in V, u is an infinite sum of generalised Fourier coefficients a_i and
x_i. And infact these generalised Fourier coefficients can be written in terms of the bilinear form c 
so that each a_i is equal to c(u, x_i).

At this point our eigenvalue problem has infinitely many solutions for lambda. And now we can look at our
oridinary differential equation. The authors show that Tn(t) has the following form given in equation 7.

Combining these solutions, we obtain a series solution for our general vibration problem. This series
solution is a linear combination of the eiegnvalues and eigenfunctions. And this shows that if we can show 
that the eigenvalues and eigenfunctions of two models are close, then the solutions of the models should 
also be close.

Also since Fourier coefficients quickly gets very small, we also only need to look at the first few eigenvalues
and eigenvectors of the models.

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

#########################################################################################################

