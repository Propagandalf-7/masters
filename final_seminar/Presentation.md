

#########################################################################################################

Thank you everyone for joining my master's seminar today. This presentation is titled 
`Finite Element Analysis of multi-dimensional and simplified models for beams`. My thesis also covers 
plate models; however, today's presentation will focus solely on beams, despite the mention of plates 
in the seminar invitation, which was an error on my part.

Beams are fundamental structural elements in engineering. They are used in a vast array of applications 
including bridges, buildings, and even the International Space Station.

Simplified models are often used over more realistic models for vibration problems due to their 
reduced complexity and lower computational demand. More realistic models are also often 
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
    to a comparison of a two and three-dimensional model. This is the main article that this presentation is \
    based on.

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

Consider the article `Comparison of linear beam theories` by Labuschagne, Van Rensburg and Van Der Merwe.

In this article, the authors investigate the efficacy of the Timoshenko beam theory against a more
realistic two-dimensional model. Specifically, the authors compare a cantilever Timoshenko beam to 
a cantilever two-dimensional beam, as a cantilever is use often in practical applications. 

The authors use modal analysis and the Finite Element Method to calculate the eigenvalues and
eigenfunctions of the models and compare these. We'll go into more detail into how the eigenvalues are
compared later in this presentation.

The authors show that the Timoshenko beam compares very well to the two-dimensional model. Especially for
the first ten eigenvalues. The authors also show that the shape of the models significantly influences 
the comparison. We'll also see this when we get to the numerical results.

#########################################################################################################

The authors do not include damping in their investigation, but show that for practical applications, the
Timshenko beam model compares well to a two-dimensional model. Practical applications do refer to common
beam shapes. 

Although a direct comparison to a three-dimensional model is preferred, the authors recommend that the
two-dimensional model be used as an intermediate step, as this should mitigate some complexities that we'll
discuss in the numerical results.

For this literature study of the article we will dicuss and replicate the results of the article, up to 
5 significant digits. We also extend the results to a comparison of a two-dimensional cantilever beam to 
a three-dimensional cantilever beam.

#########################################################################################################

Next we look at an article we'll just quickly discuss that I included for interesting results. The article
is titled `On the valid frequency range of Timoshenko beam theory`. In this article, the authors report 
on an experiment that they conducted. They suspended a small beam, made from an aluminum alloy, at both 
endpoints on carbon fibre strings.

The beam is vibrated and a free-free beam is simulated. They then sweep through the frequencies, and 
measured the natural frequencies of the beam.

The results obtained is then compared to some theoretical results. 

#########################################################################################################

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
a cantilever configuration with a retangular cross-secion of the models, following the first article 
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
For this special case, we assume that all the stress tensor components in the x3 direction is 0, the derivatives
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

Now we can start to look at the numerical results. Following what we now know about modal analysis, we 
need to solve the eigenvalue problems for our models. First we look at the eigenvalue problem for the 
Timoshenko cantilever beam.

Here we have a general eigenvalue problem for the Timoshenko beam model. Find eigenfunctions u and phi
and a real numbda lambda, the eigenvalue, satisfying equations 9 and 10. This general problem as well 
as a method to obtain a solution is presented in the article 
`Natural frequencies and modes of a Timoshenko beam` by Van Rensburg and Van Der Merwe.

In this article, the authors provide a method that allows us to calcualte the exact eigenvalues and 
eigenfunctions for the Timoshenko eigenvalue value problem. The authors also use a cantilever beam as an
example to explain the method. And then we can estimate the eigenvalues and eigenfunctions numerically
to any desired degree of accuracy.

We will not be going through the method in the presentation, and skip towards the results.


#########################################################################################################

Solving the eigenvalue problems for the two-dimensional and three-dimensional cantilever beams require a
different approach. This is where the Finite Element Method is useful.

It is not trivial to solve the eigenvalue problem for the cantilever two- and three-dimensional beams. 
The Finite Element Method (FEM) can be used to approximate the eigenvalues and eigenvectors.

We descrise our beams into a finite number of elements and nodes. For the two-dimensional beam, we can 
use a grid of rectangular shaped elements. For the three-dimensional beam, we can use a grid of cuboid 
or brick shaped elements.

We can also choose our basis functions. In this presentation, bi-cubic and tri-cubic basis function 
are used. These have the advantage over bi-linear and tri-linear basis function in that they provide
faster convergence of the eigenvalues, although they are more complex to program.


#########################################################################################################

Here we have a general form of the eigenvalue problem that is applicable to both the two and 
three-dimensional beam models that can be obtained when applying the FEM. The matrices M and K are the 
mass and stiffness matrices and are derived from the standard finite element matrices. 

In this form, the eigenvalue problem is essentially a large system of ordinary differential equations 
that can be solved using standard numerical methods.

In this presentation, we will not go through the derivation of this eigenvalue problem or the application
of the Finite Element Method, and skip towards the results.

#########################################################################################################

Since we are using a numerical method, we first look at the rate of convergence of the methods.

Here we have two figures that shows the rate of convergence of the first 20 eigenvalues of the 
two-dimensional beam on top and the three-dimensional beam at the bottom. On the ys axis is show the 
absolute differences between the eigenvalues, and on the x axis we have the number of elements needed.
Also the colors are chosen so that each color represents the same eigenvalue.

Here we see that the two-dimensional model starts to level off. For the two-dimensional model, we can get
the accuracy of the first 10 eigenvalues down to at least 5 significant digits.

The three-dimensional model actually have more room, but as you can see from the number of elements, we
started to reach a limit on what a home computer with a realistic amount of ram will be able to do. For 
the three-dimensional model we can get the first 10 eigenvalues accurate to at least 4 significant digits.

Although not shown in the graph, the shape of the beam does affect the rate of convergence.

#########################################################################################################

Now stepping back to what we want to show in this presentation. First we want to look at the validity 
of the Timoshenko beam model. Is this the work that is done in the article Comparison of linear beam 
theories. 

Here is a representation of the two beams that we are going to compare.

For the Timoshenko beam model, we do not have a paramter h representing the height of the beam. But 
something that was not shown previously in the derivation of the models, was derivation of the 
dimensionless parameter alpha. Alpha is equal to A ell squared over I where A is the area 
of the cross section, ell is the lenght of the beam and I is the area momement of inertia. And since we
have a beam with a rectangular cross-section, we can easily calculate the area and area moment of inertia.
We then obtain the following relationship between alpha and h. So for any value of h, we can use a 
corresponding value of alpha so that the beams have the same shape.

So simplicity, I will only refer to h, and the height to length ratio.

#########################################################################################################

To start, lets look at the mode shapes of the beam models. Here we have two examples of the mode shapes.
The two-dimensional cantilever beam is on the left and the Timoshenko cantilever beam is on the right.

In this example, both beams have a lenght the width ratio of 1 to 20. And we can see that the mode shapes
are similar in shape. I also included the eigenvalues, that we can also see are also quite close.

But the numbering of the eigenvalues are not the same.

#########################################################################################################

Now if we actually look at the all of the mode shapes of the two-dimensional cantilever beam and the 
Timoshenko cantilever beam, we'll see that there are some mode shapes that are present in the 
two-dimensional model that are absent in the Timoshenko beam model. 

This is due to the complexity of the model, that also introduces unqiue behaviour not found in the 
Timoshenko beam models. Now since we are interested mainly in the beam type problems, we can call the 
eigenvalies relating to the mode shapes that are shared between the models, beam-type eigenvalues and the 
rest non-beam type eigenvalues.

And we will use the mode shapes to match up our eigenvalues to be able to make our comparison.

#########################################################################################################

We can also look at the mode shapes by overlaying the two shapes. Here we have a figure that shows a mode
shapes of both the models. These shapes where only scaled to rougthly match. The scaling is not important
as any multiple of a eigenfunction is still an eigenfunction, but we are interested mainly in the overall
shape.

We can do a similar figure for the mode shape of the rotation phi. The two dimensional model does not have
mode shapes representing the angle of the cross-sections. But be can estimate the rotation by calculating
the average rotation of the cross-sections at intervals. This is demonstrated in Figure 6 by the red line.
And we actually get quite a good approximation.


#########################################################################################################

Now lets look at some results from the article Comparison of linear beam models.

On the left of the table, we have the results from the article and on the right we have the same results,
replicated in my master dissertation.

These results show that the eigenvalues match quite good, with a maximum relative error less than 3%. 
These results are for a beam with height to lenght ratio of 1 to 10. 

#########################################################################################################

Lets look at similar results for different shapes of beams. 

In this table, we have the first 20 eigenvalues of the two-dimensional beam model matched to the eiegnvalues
of the Timoshenko beam model for different beam shapes. 

From left to right the height to length ratio decreases as the beam goes from shot and thick to long 
and slender. And in these results, we see that the maximum relative error is decreasing as the beam 
gets more slender. This result is also show in the article by hands of an illustration.

Interestingly we also see that the number of non-beam type eigenvalues decrease as the beam gets more 
slender as we could match 13 of the eigenvalues for the a height to length ratio of 1 to 5, and 15 eigenvalues
for the a height to length ratio of 1 to 30.


Overall we see that the Timoshenko beam model compares very well to the Two dimensional beam model. The 
shape of the beam influences how well the model compares. We were also able to confirm the results of 
the article.

Next we extend this to investigate the validity of the two-dimensional model.


#########################################################################################################

Next we look at the comparison of the two-dimensional cantilever beam and a three-dimensional catilever 
beam. This is an extention of the article.

Again we start be looking at the mode shapes. Here we have two examples of mode shapes relating to beam
type eigenvalues. Again the shapes are very similar as are the eigenvalues. We also again see a difference 
between the numbering of the eigenvalues, again indicating a difference in the complexity of the models.

And it is exactly due to this complexity that the authors suggest using the two-dimensional model as 
an intermeidate step and not directly comparing the one-dimensional Timoshenko model to the
three-dimensional model.

#########################################################################################################

Here we have some of the non-beam type mode shapes. The first row shows a mode shape that is shared between
the two models. In fact all of the mode shapes of the two-dimensional model are present in the 
three-dimensional model, but this is expected.

The second row shows examples of mode shapes that are unique to the three-dimensional model. These mode
shapes represent interesting behaviour that the three-dimensional model can extert. 

In fact there are quite a lot of extras, so in the results that follow, we will only include the 
eigenvalues that are shared between the models, with the non-beam type eigenvalues grayed out.

#########################################################################################################

We also have an extra parameter for the depth of the beam, so we will start by looking at a short and 
thick beam with a lenght to width ration of 1 to 5, and then vary the depth of the beam, first looking
at a depth that is smaller than the hieight of the beam.

In this table we have the results of the two-dimensional beam on the far right hand side. and then from 
left to right, the depth of the beam is decreased. Remeber that the non-beam type eigenvalues are 
grayed out.

Here we can see that as the depth of the beam decreases, our maximum relative error also decreases. 
Although overall our comparison is very good, with all results less that a 2.7% difference. These 
maximum relative errors do include the non-beam type results, but if we split them,


#########################################################################################################

We get the following table. Here we just show the maximum relative errors distinctly for the beam type 
and non beam type eigenvalues of the previous table. For beam type problems, the two-dimensional model is
doing very well.

#########################################################################################################

In this table, we have similar results as the previous one, but now for a long and slender beam. Recall 
that the long and slender beam was the best case in the comparison of the Timoshenko cantilever beam and 
the two-dimensional cantilever beam.

Here again we see that the long and slender beam compares better than the shorter and thick beam. And
as we decrease the depth of the beam, we get exceptional results.

Something to note that I did not mention when we looked at the previous table, is that as the depth of 
the beam decreases, the number of non-beam type eigenvalues increases. This is different to what we see
when we decrease the height to lenght ratio.

This also shows how computationally expensive the three-dimensional model is and why a simpler model
is often used, especially for beam type problems.

#########################################################################################################

Here again we split the beam type and non beam types results.  

#########################################################################################################

Next we look at what happens if we increase the depth of the beam. Here we have a table showing a beam 
of length to width ratio of 1 to 20 (which was the best case in all of our results), and we increase the 
depth of the beam.

We do not need a lot of results to see that as the depth increases, our maximum relative error increases.
We also again see an increase of non-beam type eigenvalues, suggesting any variation from a square beam 
will be a little bit more complex.

#########################################################################################################

Again we split the results up into beam type and non-beam type and we wee that as we increase the 
depth, our beam type comparison get worse.

Now this is not necessarily a fault with the two-dimensional beam model. It suggests that for 
this shapes, a beam model is not appriopriate. A model like a plate model might be better suited. In my 
masters dissertation, I did expand the investigation to plate models and it shows that a plate model
is well suited if the depth is larger than the height of the body.

#########################################################################################################

#########################################################################################################

#########################################################################################################