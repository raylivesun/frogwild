# FrogWild! – Fast PageRank Approximations
# on Graph Engines
    
<html xmlns="http://www.w3.org/1999/xhtml" lang="en"><head>
<meta http-equiv="content-type" content="text/html; charset=UTF-8">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
       </head>
    <body>
        <div>
            <h1>FrogWild! – Fast PageRank Approximations
                on Graph Engines</h1>    
        </div>
        <div align="center">
            <p>Ioannis Mitliagkas <br>
                Michael Borokhovich <br>
                ECE, UT AustinECE, UT Austin <br>
                ioannis@utexas.edu <br>
                Alexandros G. Dimakismichaelbor@utexas.edu <br>
                Constantine Caramanis <br>
                ECE, UT AustinECE, UT Austin <br>
                dimakis@austin.utexas.edu constantine@utexas.edu <br>
                ABSTRACT<br>
            </p>                 
        </div>
        <div>
            <table width="100%" height="100%">
                <tbody><tr>
                    <td>
                        <h2>ABSTRACT</h2>                       
                         <p>We propose FrogWild, a novel algorithm for fast approxi-
                            mation of high PageRank vertices, geared towards reducing
                            network costs of running traditional PageRank algorithms.
                            Our algorithm can be seen as a quantized version of power
                            iteration that performs multiple parallel random walks over
                            a directed graph. One important innovation is that we in-
                            troduce a modification to the GraphLab framework that
                            only partially synchronizes mirror vertices. This partial
                            synchronization vastly reduces the network traffic generated
                            by traditional PageRank algorithms, thus greatly reducing
                            the per-iteration cost of PageRank. On the other hand,
                            this partial synchronization also creates dependencies be-
                            tween the random walks used to estimate PageRank. Our
                            main theoretical innovation is the analysis of the correla-
                            tions introduced by this partial synchronization process and
                            a bound establishing that our approximation is close to the
                            true PageRank vector.</p>
                            <h2>1.
                                INTRODUCTION</h2>
                                <p>Large-scale graph processing is becoming increasingly im-
                                    portant for the analysis of data from social networks, web
                                    pages, bioinformatics and recommendation systems. Graph
                                    algorithms are difficult to implement in distributed com-
                                    putation frameworks like Hadoop MapReduce and Spark.
                                    For this reason several in-memory graph engines like Pregel,
                                    Giraph, GraphLab and GraphX [24, 23, 35, 31] are being
                                    developed. There is no full consensus on the fundamental
                                    abstractions of graph processing frameworks but certain pat-
                                    terns such as vertex programming and the Bulk Synchronous
                                    Parallel (BSP) framework seem to be increasingly popular.
                                    PageRank computation [27], which gives an estimate of
                                    the importance of each vertex in the graph, is a core compo-
                                    nent of many search routines; more generally, it represents,
                                    de facto, one of the canonical tasks performed using such
                                    graph processing frameworks. Indeed, while important in
                                    its own right, it also represents the memory, computation
                                    and communication challenges to be overcome in large scale
                                    iterative graph algorithms.</p>
                                    <p>user activity graph), PageRank should be recalculated con-
                                        stantly. Moreover, the key users constitute only a small
                                        fraction of the total number of users, thus, a fast approxi-
                                        mation for the top-PageRank nodes constitutes a desirable
                                        alternative to the exact solution.
                                        In this paper we address this problem. Our algorithm
                                        (called FrogWild for reasons that will become subsequently
                                        apparent) significantly outperforms the simple reduced iter-
                                        ations heuristic in terms of running time, network communi-
                                        cation and scalability. We note that, naturally, we compare
                                        our algorithm and reduced-iteration-PageRank within the
                                        same framework: we implemented our algorithm in GraphLab
                                        PowerGraph and compare it against the built-in PageRank
                                        implementation. A key part of our contribution also involves
                                        the proposal of what appears to be simply a technically mi-
                                        nor modification within the GraphLab framework, but nev-
                                        ertheless results in significant network-traffic savings, and
                                        we believe may nevertheless be of more general interest be-
                                        yond PageRank computations.</p>
                                        <p>Contributions: We consider the problem of fast and
                                            efficient (in the sense of time, computation and communica-
                                            tion costs) computation of the high PageRank nodes, using
                                            a graph engine. To accomplish this we propose and ana-
                                            lyze an new PageRank algorithm specifically designed for
                                            the graph engine framework, and, significantly, we propose
                                            a modification of the standard primitives of the graph en-
                                            gine framework (specifically, GraphLab PowerGraph), that
                                            enables significant network savings. We explain in further
                                            detail both our objectives, and our key innovations.
                                            Rather than seek to recover the full PageRank vector, we
                                            aim for the top k PageRank vertices (where k is considered
                                            to be approximately in the order of 10 − 1000). Given an
                                            output of a list of k vertices, we define two natural accuracy
                                            metrics that compare the true top-k list with our output.
                                            The algorithm we propose, FrogWild operates by start-
                                            ing a small (sublinear in the number of vertices n) number
                                            of random walkers (frogs) that jump randomly on the di-
                                            rected graph. The random walk interpretation of PageRank
                                            enables the frogs to jump to a completely random vertex
                                            (teleport) with some constant probability (set to 0.15 in our
                                            experiments, following standard convention). After we al-
                                            low the frogs to jump for time equal to the mixing time of
                                            this non-reversible Markov chain, their positions are sam-
                                            pled from the invariant distribution π which is normalized
                                            PageRank. The standard PageRank iteration can be seen as
                                            the continuous limit of this process (i.e., the frogs become
                                            water), which is equivalent to power iteration for stochastic
                                            matrices.</p>
                                            <p>standard analysis of the Power Method iteration no longer
                                                applies. The main challenge that arises is the theoretical
                                                analysis of the FrogWild algorithm. The model is that
                                                each vertex is separated across machines and each connec-
                                                tion between two vertex copies is present with probability ps .
                                                A single frog performing a random walk on this new graph
                                                defines a new Markov Chain and this can be easily designed
                                                to have the same invariant distribution π equal to normal-
                                                ized PageRank. The complication is that the trajectories
                                                of frogs are no longer independent: if two frogs are in ver-
                                                tex v1 and (say) only one mirror v10 synchronizes, both frogs
                                                will need to jump through edges connected with that par-
                                                ticular mirror. Worse still, this correlation effect increases,
                                                the more we seek to improve network traffic by further de-
                                                creasing ps . Therefore, it is no longer true that one obtains
                                                independent samples from the invariant distribution π. Our
                                                theoretical contribution is the development of an analytical
                                                bound that shows that these dependent random walks still
                                                can be used to obtain π̂ that is provably close to π with high
                                                probability. We rely on a coupling argument combined with
                                                an analysis of pairwise intersection probabilities for random
                                                walks on graphs. In our convergence analysis we use the
                                                contrast bound [12] for non-reversible chains.</p>
                                                <p>right eigenvector of Q. That is, π , v1 (Q). By the Perron-
Frobenius theorem, the corresponding eigenvalue is 1. This
implies the fixed-point characterization of the PageRank vec-
tor, π = Qπ.
The PageRank vector assigns high values to important
nodes. Intuitively, important nodes have many important
predecessors (other nodes that point to them). This recur-
sive definition is what makes PageRank robust to manipula-
tion, but also expensive to compute. It can be recovered by
exact eigendecomposition of Q, but at real problem scales
this is prohibitively expensive. In practice, engineers often
use a few iterations of the power method to get a ”good-
enough” approximation.
The definition of PageRank hinges on the left-stochastic
matrix Q, suggesting a connection to Markov chains. In-
deed, this connection is well documented and studied [2,
16]. An important property of PageRank from its random
walk characterization, is the fact that π is the invariant dis-
tribution for a Markov chain with dynamics described by
Q. A non-zero pT , also called the teleportation probability,
introduces a uniform component to the PageRank vector π.
We see in our analysis that this implies ergodicity and faster
mixing for the random walk.</p>
<h2>1.1 Notation</h2>
<p>Lowercase letters denote scalars or vectors. Uppercase
    letters denote matrices. The (i, j) element of a matrix A
    is Aij . We denote the transpose of a matrix A by A0 . For
    a time-varying vector x, we denote its value at time t by
    xt . When not otherwise specified, kxk denotes the l2 -norm
    of vector x. We use ∆n−1 for the probability simplex in
    n dimensions, and and ei ∈ ∆n−1 for the indicator vector
    for item i. For example, e1 = [1, 0, ...0]. For the set of all
    integers from 1 to n we write [n].</p>
    <p>2.2
        Algorithm</p>
        <p>During setup, the graph is partitioned using GraphLab’s
            default ingress algorithm. At this point each one of N frogs
            is born on a vertex chosen uniformly at random. Each vertex
            i carries a counter initially set to 0 and denoted by c(i).
            Scheduled vertices execute the following program.
            Incoming frogs from previously executed vertex programs,
            are collected by the init() function. At apply() every frog
            dies with probability pT = 0.15. This, along with a uniform
            starting position, effectively simulates the 15% uniform com-
            ponent from Definition 1.
            A crucial part of our algorithm is the change in synchro-
            nization behaviour. The <sync> step only synchronizes a
            ps fraction of mirrors leading to commensurate gains in net-
            work traffic (cf. Section 3). This patch on the GraphLab
            codebase was only a few lines of code. Section 3 contains
            more details regarding the implementation.
            The scatter() phase is only executed for edges e incident
            to a mirror of i that has been synchronized. Those edges
            draw a binomial number of frogs to send to their other end-
            point. The rest of the edges perform no computation. The
            frogs sent to vertex j at the last step will be collected at the
            init() step when j executes.</sync></p>
            <p>Here c(i) refers to the tally maintained by the FrogWild!
                vertex program.
                Assuming θ = 2.2 and picking, for example, γ = 0.5, we get
                that
                √
                P(kπk∞ &gt; 1/ n) ≤ cn−1/3 .
                Now we can state the main result. Here we give a guaran-
                tee for the quality of the solution furnished by our algorithm.
                This implies that with probability at least 1 − cn−1/3 the
                meeting probability is bounded as follows.
                Theorem 1 (Main Theorem). Consider N frogs fol-
                lowing the FrogWild! process (Section 2.2), under the era-
                sure model of Definition 8. The frogs start at independent
                locations, distributed uniformly and stop after a geometric
                number of steps or, at most, t steps. The estimator π̂N
                (Definition 5), captures mass close to the optimal. Specifi-
                cally, with probability at least 1 − δ,
                p∩ (t) ≤
                One would usually take a number of steps t that are either
                constant or logarithmic with respect to the graph size n.
                This implies that for many reasonable choices of set size k
                and acceptable probability of failure δ, the meeting proba-
                bility vanishes as n grows. Then we can make the second
                term of the error in (4) arbitrarily small by controlling the
                number of frogs, N . The proof for Proposition 7 is deferred
                to Appendix B.3.
                µk (π̂N ) ≥ µk (π) −,</p>
                <p>Remark 6 (Scaling). The result in Theorem 1 imme-
                    diately implies the following scaling for the number of itera-
                    tions and frogs respectively. They both depend on the max-
                    imum captured mass possible, µk (π) and are sufficient for
                    making the error, , of the same order as µk (π).</p>
                    <p>The proof of Theorem 1 is deferred to Appendix B.1.
                        The guaranteed accuracy via this result also depends on
                        the probability that two walkers will intersect. Via a simple
                        argument, that probability is the same as the meeting prob-
                        ability for independent walks. The next theorem calculates
                        this probability.
                        Theorem 2 (Intersection Probability). Consider
                        two independent random walks obeying the same ergodic tran-
                        sition probability matrix, Q with invariant distribution π, as
                        described in Definition 1. Furthermore, assume that both of
                        them are initially distributed uniformly over the state space
                        of size n. The probability that they meet within t steps, is
                        bounded as follows,</p>
                        <p>where kπk∞ , denotes the maximal element of the vector π.
                            The proof is based on the observation that the l∞ norm of
                            a distribution controls the probability that two independent
                            samples coincide. We show that for all steps of the random
                            walk, that norm is controlled by the l∞ norm of π. We defer
                            the full proof to Appendix B.2.
                            A number of studies, give experimental evidence (e.g. [8])
                            suggesting that PageRank values for the web graph follow
                            a power-law distribution with parameter approximately θ =
                            2.2. That is true for the tail of the distribution – the largest
                            values, hence of interest to us here – regardless of the choice
                            of pT . The following proposition bounds the value of the
                            heaviest PageRank value, kπk∞ .
                            Proposition 7 (Max of Power-Law Distribution).
                            Let π ∈ ∆n−1 follow a power-law distribution with param-
                            eter θ and minimum value pT /n. Its maximum element,
                            1
                            kπk∞ , is at most n−γ , with probability at least 1 − cnγ− θ−1 ,
                            for some universal constant c.</p>
                    </td><td>
                        <p>In this paper we propose a novel algorithm for fast approx-
                            imate calculation of high PageRank vertices. Note that even
                            though most previous works calculate the complete PageR-
                            ank vector (of length in the millions or billions), in many
                            graph analytics scenarios a user wants a quick estimation of
                            the most important or relevant nodes – distinguishing the
                            10th most relevant node from the 1, 000th most relevant is
                            important; the 1, 000, 000th from the 1, 001, 000th much less
                            so. A simple solution is to run the standard PageRank algo-
                            rithm for fewer iterations (or with an increased tolerance).
                            While certainly incurring less overall cost, the per-iteration
                            cost remains the same; more generally, the question remains
                            whether there is a more efficient way to approximately re-
                            cover the heaviest PageRank vertices.
                            There are many real life applications that may benefit
                            from a fast top-k PageRank algorithm. One example is grow-
                            ing loyalty of influential customers [1]. In this application,
                            a telecom company identifies the top-k influential customers
                            using the top-k PageRank on the customers’ activity (e.g.,
                            calls) graph. Then, the company invests its limited bud-
                            get on improving user experience for these top-k customers,
                            since they are most important for building good reputation.
                            Another interesting example is an application of PageRank
                            for finding keywords and key sentences in a given text. In
                            [25], the authors show that PageRank performs better than
                            known machine learning techniques for keyword extraction.
                            Each unique word (noun, verb or an adjective) is regarded
                            as a vertex and there is an edge between two words if they
                            occur in close proximity in the text. Using approximate top-
                            k PageRank, we can identify the top-k keywords much faster
                            than obtaining the full ranking. When keyword extraction is
                            used by time sensitive applications or for an ongoing analysis
                            of a large number of documents, speed becomes a crucial fac-
                            tor. The last example we describe here is the application of
                            PageRank for online social networks (OSN). It is important
                            in the context of OSNs to be able to predict which users will
                            remain active in the network for a long time. Such key users
                            play a decisive role in developing effective advertising strate-
                            gies and sophisticated customer loyalty programs, both vital
                            for generating revenue [19]. Moreover, the remaining users
                            can be leveraged, for instance for targeted marketing or pre-
                            mium services. It is shown in [19] that PageRank is a much
                            more efficient predictive measure than other centrality mea-
                            sures. The main innovation of [19] is the usage of a mixture
                            of connectivity and activity graphs for PageRank calcula-
                            tion. Since these graphs are highly dynamic (especially the</p>
                            <p>One major issue with simulating discrete frogs in a graph
                                engine is teleportations. Graph frameworks partition ver-
                                tices to physical nodes and restrict communication on the
                                edges of the underlying graph. Global random jumps would
                                create dense messaging patterns that would increase com-
                                munication. Our second innovation is a way of obtaining
                                an identical sampling behavior without teleportations. We
                                achieve this by initiating the frogs at uniformly random posi-
                                tions and having them perform random walks for a life span
                                that follows a geometric random variable. The geometric
                                probability distribution depends on the teleportation prob-
                                ability and can be calculated explicitly.
                                Our third innovation involves a simple proposed modifica-
                                tion for graph frameworks. Most modern graph engines (like
                                GraphLab PowerGraph [17]) employ vertex-cuts as opposed
                                to edge-cuts. This means that each vertex of the graph is
                                assigned to multiple machines so that graph edges see a local
                                vertex mirror. One copy is assigned to be the master and
                                maintains the master version of vertex data while remaining
                                replicas are mirrors that maintain local cached read–only
                                copies of the data. Changes to the vertex data are made to
                                the master and then replicated to all mirrors at the next syn-
                                chronization barrier. This architecture is highly suitable for
                                graphs with high-degree vertices (as most real-world graphs
                                are) but has one limitation when used for a few random
                                walks: imagine that vertex v1 contains one frog that wants
                                to jump to v2 . If vertex v1 has very high degree, it is very
                                likely that multiple replicas of that vertex exist, possibly
                                one in each machine in the cluster. In an edge-cut scenario
                                only one message would travel from v1 → v2 , assuming v1
                                and v2 are located in different physical nodes. However,
                                when vertex-cuts are used, the state of v1 is updated (i.e.,
                                contains no frogs now) and this needs to be communicated
                                to all mirrors. It is therefore possible that a single random
                                walk can create a number of messages equal to the number
                                of machines in the cluster.
                                We modify PowerGraph to expose a scalar parameter ps
                                per vertex. By default, when the framework is running, in
                                each super-step all masters synchronize their programs and
                                vertex data with their mirrors. Our modification is that
                                for each mirror we flip an independent coin and synchronize
                                with probability ps . Note that when the master does not
                                synchronize the vertex program with a replica, that replica
                                will not be active during that super-step. Therefore, we can
                                avoid the communication and CPU execution by performing
                                limited synchronization in a randomized way.
                                FrogWild is therefore executed asynchronously but re-
                                lies on the Bulk Synchronous execution mode of PowerGraph
                                with the additional simple randomization we explained. The
                                name of our algorithm is inspired by HogWild [29], a lock-
                                free asynchronous stochastic gradient descent algorithm pro-
                                posed by Niu et al.. We note that PowerGraph does support
                                an asynchronous execution mode [17] but we implemented
                                our algorithm by a small modification of synchronous execu-
                                tion. As discussed in [17], the design of asynchronous graph
                                algorithms is highly nontrivial and involves locking proto-
                                cols and other complications. Our suggestion is that for the
                                specific problem of simulating multiple random walks on a
                                graph, simply randomizing synchronization can give signifi-
                                cant benefits while keeping design simple.
                                While the parameter ps clearly has the power to signifi-
                                cantly reduce network traffic – and indeed, this is precisely
                                born out by our empirical results – it comes at a cost: the</p>
                                <p>right eigenvector of Q. That is, π , v1 (Q). By the Perron-
                                    Frobenius theorem, the corresponding eigenvalue is 1. This
                                    implies the fixed-point characterization of the PageRank vec-
                                    tor, π = Qπ.
                                    The PageRank vector assigns high values to important
                                    nodes. Intuitively, important nodes have many important
                                    predecessors (other nodes that point to them). This recur-
                                    sive definition is what makes PageRank robust to manipula-
                                    tion, but also expensive to compute. It can be recovered by
                                    exact eigendecomposition of Q, but at real problem scales
                                    this is prohibitively expensive. In practice, engineers often
                                    use a few iterations of the power method to get a ”good-
                                    enough” approximation.
                                    The definition of PageRank hinges on the left-stochastic
                                    matrix Q, suggesting a connection to Markov chains. In-
                                    deed, this connection is well documented and studied [2,
                                    16]. An important property of PageRank from its random
                                    walk characterization, is the fact that π is the invariant dis-
                                    tribution for a Markov chain with dynamics described by
                                    Q. A non-zero pT , also called the teleportation probability,
                                    introduces a uniform component to the PageRank vector π.
                                    We see in our analysis that this implies ergodicity and faster
                                    mixing for the random walk.</p>
                                    <p>right eigenvector of Q. That is, π , v1 (Q). By the Perron-
                                        Frobenius theorem, the corresponding eigenvalue is 1. This
                                        implies the fixed-point characterization of the PageRank vec-
                                        tor, π = Qπ.
                                        The PageRank vector assigns high values to important
                                        nodes. Intuitively, important nodes have many important
                                        predecessors (other nodes that point to them). This recur-
                                        sive definition is what makes PageRank robust to manipula-
                                        tion, but also expensive to compute. It can be recovered by
                                        exact eigendecomposition of Q, but at real problem scales
                                        this is prohibitively expensive. In practice, engineers often
                                        use a few iterations of the power method to get a ”good-
                                        enough” approximation.
                                        The definition of PageRank hinges on the left-stochastic
                                        matrix Q, suggesting a connection to Markov chains. In-
                                        deed, this connection is well documented and studied [2,
                                        16]. An important property of PageRank from its random
                                        walk characterization, is the fact that π is the invariant dis-
                                        tribution for a Markov chain with dynamics described by
                                        Q. A non-zero pT , also called the teleportation probability,
                                        introduces a uniform component to the PageRank vector π.
                                        We see in our analysis that this implies ergodicity and faster
                                        mixing for the random walk.</p>
                                        <h2>2.1.1 Top PageRank Element</h2>
                                        <p>Given the true PageRank vector, π and an estimate v
                                            given by an approximate PageRank algorithm, we define
                                            the top-k accuracy using two metrics.</p>
                                            <p>Definition 2 (Mass Captured). Given a distribution
                                                v ∈ ∆n−1 , the true PageRank distribution π ∈ ∆n−1 and an
                                                integer k ≥ 0, we define the mass captured by v as follows.</p>
                                                <p>Put simply, the set S ∗ that gets the most mass according
                                                    to v out of all sets of size k, is evaluated according to π and
                                                    that gives us our metric. It is maximized by π itself, i.e. the
                                                    optimal value is µk (π).
                                                    The second metric we use is the exact identification prob-
                                                    ability, i.e. the fraction of elements in the output list that
                                                    are also in the true top-k list. Note that the second metric is
                                                    limited in that it does not give partial credit for high PageR-
                                                    ank vertices that were not in the top-k list. In our experi-
                                                    ments in Section 3, we mostly use the normalized captured
                                                    mass accuracy metric but also report the exact identification
                                                    probability for some cases – typically the results are similar.
                                                    We subsequently describe our algorithm. We attempt to
                                                    approximate the heaviest elements of the invariant distribu-
                                                    tion of a Markov Chain, by simultaneously performing mul-
                                                    tiple random walks on the graph. The main modification to
                                                    PowerGraph, is the exposure of a parameter, ps , that con-
                                                    trols the probability that a given master node synchronizes
                                                    with any one of its mirrors. Per step, this leads to a propor-
                                                    tional reduction in network traffic. The main contribution
                                                    of this paper is to show that we get results of comparable
                                                    or improved accuracy, while maintaining this network traffic
                                                    advantage. We demonstrate this empirically in Section 3.</p>
                    <p>denoted by X, is geometrically distributed with parameter
                        pT . This follows from the time-reversibility in the telepor-
                        tation process: inter-teleportation times are geometrically
                        distributed, so as long as the first teleportation event hap-
                        pens before the stopping time, then X ∼ Geom(pT ).
                        This establishes that, the FrogWild! process – where a
                        frog performs a geometrically distributed number of steps
                        following the original transition matrix P – closely mimics
                        a random walk that follows the adjusted transition matrix,
                        Q. In practice, we stop the process after t steps to get a
                        good approximation. To show our main result, Theorem 1,
                        we analyze the latter process.
                        Using a binomial distribution to independently generate
                        the number of frogs in the scatter() phase closely mod-
                        els the effect of random walks. The marginal distributions
                        are correct, and the number of frogs, that did not die dur-
                        ing the apply() step, is preserved in expectation. For our
                        implementation we resort to a more efficient approach. As-
                        suming K(i) frogs survived the apply() step, and M mirrors
                        e frogs
                        where picked for synchronization, then we send d K(i)
                        M
                        to min(K(i), M ) mirrors. If the number of available frogs is
                        less than the number of synchronized mirrors, we pick K(i)
                        arbitrarily.</p> 
                        <hr>
                        <h2>FrogWild! vertex program</h2>
                        <h2>2.3
                            FrogWild! vertex program
                            Main Result</h2>
                            <p>Our analytical results essentially provide a high probabil-
                                ity guarantee that our algorithm produces a solution that
                                approximates well the PageRank vector. Recall that the
                                main modification of our algorithm involves randomizing the
                                synchronization between master nodes and mirrors. For our
                                analysis, we introduce a broad model to deal with partial
                                synchronization, in Appendix A.
                                Our results tell us that partial synchronization does not
                                change the distribution of a single random walk. To make
                                this and our other results clear, we need the simple defini-
                                tion.</p>
                                <p>Definition 3. We denote the state of random walk i at
                                    its tth step by sti .
                                    Then, we see that P st+1
                                    = i st1 = j = 1/dout (j), and
                                    1
                                    xt+1
                                    = P xt1 . This follows simply by the symmetry assumed
                                    1
                                    in Definition 8. Thus if we were to sample in serial, the
                                    modification of the algorithm controlling (limiting) synchro-
                                    nization would not affect each sample, and hence would not
                                    affect our estimate of the invariant distribution. However,
                                    we start multiple (all) random walks simultaneously. In this
                                    setting, the fundamental analytical challenge stems from the
                                    fact that any set of random walks with intersection are now
                                    correlated. The key to our result is that we can control the
                                    effect of this correlation, as a function the parameter ps and
                                    the pairwise probability that two random walks intersect. We
                                    define this formally.</p>
                                    <p>Definition 4. Suppose two walkers l1 and l2 start at the
                                        same time and perform t steps. The probability that they
                                        meet is defined as follows.
                                        p∩ (t) , P (∃ τ ∈ [0, t], s.t. sτl1 = sτl2 )
                                        (2)
                                        Definition 5 (Estimator). Given the positions of N
                                        random walks at time t, {stl }N
                                        l=1 , we define the following
                                        estimator for the invariant distribution π.
                                        π̂N (i) ,
                                        4
                                        {l : l ∈ [N ], stl = i}
                                        c(i)
                                        =
                                        N
                                        N</p>
                         <p>4
                            Remark 6 (Scaling). The result in Theorem 1 imme-
                            diately implies the following scaling for the number of itera-
                            tions and frogs respectively. They both depend on the max-
                            imum captured mass possible, µk (π) and are sufficient for
                            making the error, , of the same order as µk (π).
                            
                            
                            
                            
                            1
                            k
                            t = O log
                            ,
                            N =O
                            µk (π)
                            µk (π)2
                            Prior Work
                            There is a very large body of work on computing and
                            approximating PageRank on different computation models
                            (e.g. see [10, 13, 30, 14, 4] and references therein). To the
                            best of our knowledge, our work is the first to specifically de-
                            sign an approximation algorithm for high-PageRank nodes
                            for graph engines. Another line of work looks for Personal-
                            ized PageRank (PPR) scores. This quantifies the influence
                            an arbitrary node i has on another node j, cf. recent work
                            [22] and discussion therein. In [6], the top-k approxima-
                            tion of PPR is studied. However, PPR is not applicable in
                            our case, as we are looking for an answer close to a global
                            optimum.
                            In [5], a random-walks-based algorithm is proposed. The
                            authors provide some insightful analysis of different varia-
                            tions of the algorithm. They show that starting a single
                            walker from every node, is sufficient to achieve a good global
                            approximation. We focus on capturing a few nodes with a lot
                            of mass, hence we can get away with orderwise much fewer
                            frogs than O(n). This is important for achieving low net-
                            work traffic when the algorithm is executed on a distributed
                            graph framework. Figure 8 shows linear reduction in net-
                            work traffic when the number of initial walkers decreases.
                            Furthermore, our method does not require waiting for the
                            last frog to naturally expire (note that the geometric distri-
                            bution has infinite support). We impose a very short time
                            cut-off, t, and exactly analyze the penalty in captured mass
                            we pay for it in Theorem 1.
                            One natural question is how our algorithm compares to,
                            or can be complemented by, graph sparsification techniques.
                            One issue here is that graph sparsification crucially depends
                            on the similarity metric used. Well-studied properties that
                            are preserved by different sparsification methods involve lengths
                            of shortest paths between vertices (such sparsifiers are called
                            Spanners, see e.g. [28]), cuts between subsets of vertices [9]
                            and more generally quadratic forms of the graph laplacian
                            [33, 7], see [7] and references therein for a recent overview.
                            To the best of our knowledge, there are no known graph
                            sparsification techniques that preserve vertex PageRank.
                            One natural heuristic that one may consider is to inde-
                            pendently flip a coin and delete each edge of the graph with
                            some probability r. Note that this is crucially different from
                            spectral sparsifiers [33, 7] that choose these probabilities us-
                            ing a process that is already more complicated than esti-
                            mating PageRank. This simple heuristic of independently</p>
                            <p>deleting edges indeed accelerates the estimation process for
                                high-PageRank vertices. We compare FrogWild to this
                                uniform sparsification process in Figure 5. We present here
                                results for 2 iterations of the GraphLab PR on the spar-
                                sified graph. Note that running only one iteration is not
                                interesting since it actually estimates only the in-degree of
                                a node which is known in advance (i.e., just after the graph
                                loading) in a graph engine framework. It can be seen in
                                Figure 5 that even when only two iterations are used on
                                the sparsified graph the running time is significantly worse
                                compared to FrogWild and the accuracy is comparable.
                                Our base-line comparisons come from the graph frame-
                                work papers since PageRank is a standard benchmark for
                                running-time, network and other computations. Our imple-
                                mentation is on GraphLab (PowerGraph) and significantly
                                outperforms the built-in PageRank algorithm. This algo-
                                rithm is already shown in [17, 31] to be significantly more
                                efficient compared to other frameworks like Hadoop, Spark,
                                Giraph etc.</p>
                                <h2>3. Experiments</h2>
                                <p>In this section we compare the performance of our algo-
                                    rithm to the PageRank algorithm shipped with GraphLab
                                    v2.2 (PowerGraph) [23]. The fact that GraphLab is the
                                    fastest distributed engine for PageRank is established ex-
                                    perimentally in [31]. We focus on two algorithms: the basic
                                    built-in algorithm provided as part of the GraphLab graph
                                    analytics toolkit, referred to here as GraphLab PR, and
                                    FrogWild. Since we are looking for a top-k approximation
                                    and GraphLab PR is meant to find the entire PageRank
                                    vector, we only run it for a small number of iterations (usu-
                                    ally 2 are sufficient). This gives us a good top-k approxi-
                                    mation and is much faster than running the algorithm until
                                    convergence. We also fine tune the algorithm’s tolerance
                                    parameter to get a good but fast approximation.
                                    We compare several performance metrics, namely: run-
                                    ning time, network usage, and accuracy. The metrics do
                                    not include time and network usage required for loading the
                                    graph into GraphLab (known as the ingress time). They
                                    reflect only the execution stage.</p>
                                    <h2>3.1
                                        The Systems</h2>
                                        <p>We perform experiments on two systems. The first system
                                            is a cluster of 20 virtual machines, created using VirtualBox
                                            4.3 [34] on a single physical server. The server is based on an
                                            Intel R Xeon R CPU E5-1620 with 4 cores at 3.6 GHz, and
                                            16 GB of RAM. The second system, comprises of a cluster
                                            of up to 24 EC2 machines on AWS (Amazon web services)
                                            [3]. We use m3.xlarge instances, based on Intel R Xeon R
                                            CPU E5-2670 with 4 vCPUs and 15 GB RAM.</p>
                                            <h2>3.2
                                                The Data</h2>
                                                <p>For the VirtualBox system, we use the LiveJournal graph
                                                    [21] with 4.8M vertices and 69M edges. For the AWS system,
                                                    in addition to the LiveJournal graph, we use the Twitter
                                                    graph [20] which has 41.6M nodes and 1.4B edges.</p>
                                                    <p>3.3 implementation</p>
                                                    <p>FrogWild is implemented on the standard GAS (gather,
                                                        apply, scatter) model. We implement init(), apply(), and
                                                        scatter(). The purpose of init() is to collect the random
                                                        walks sent to the node by its neighbors using scatter() in</p>
                     <p>he previous iteration. In the first iteration, init() gener-
                        ates a random fraction of the initial total number of walkers.
                        This implies that the initial walker locations are randomly
                        distributed across nodes. FrogWild requires the length of
                        random walks to be geometrically distributed (see Section
                        2.2). For the sake of efficiency, we impose an upper bound on
                        the length of random walks. The algorithm is executed for
                        the constant number of iterations (experiments show good
                        results with even 3 iterations) after which all the random
                        walks are stopped simultaneously. The apply() function
                        is responsible for keeping track of the number of walkers
                        that have stopped on each vertex and scatter() distributes
                        the walkers still alive to the neighbors of the vertex. The
                        scatter() phase is the most challenging part of the imple-
                        mentation. In order to reduce information exchange between
                        machines, we use a couple of ideas.
                        First, we notice that random walks do not have iden-
                        tity. Hence, random walks destined for the same neighbor
                        can be combined into a single message. The second opti-
                        mization and significant part of our work is modifying the
                        GraphLab engine. The recent versions of GraphLab (since
                        PowerGraph) partition the graph by splitting vertices. As a
                        consequence, the engine will need to synchronize all the mir-
                        rors of a vertex over the network a number of times during
                        each GAS cycle.
                        When running a few random walks, only a handful of
                        neighbors end up receiving walkers. For this reason, syn-
                        chronizing all mirrors can be very wasteful. We deal with
                        that by implementing randomized synchronization. We ex-
                        pose parameter ps ∈ [0, 1] to the user as a small extension to
                        the GraphLab API. It describes the fraction of replicas that
                        will be synchronized. Replicas not synchronized remain idle
                        for the upcoming scatter phase. The above edits in the en-
                        gine are only a matter of a few (about 10) lines of code. Note
                        that the ps parameter is completely optional, i.e., setting it
                        to 1 will result in the original engine operation. Hence, other
                        analytic workloads will not be affected. However, any ran-
                        dom walk or “gossip” style algorithm (that sends a single
                        messages to a random subset of its neighbors) can benefit
                        by exploiting ps . Our modification of the GraphLab engine
                        as well as the FrogWild vertex program can be found in
                        [11].</p>
                        <h2>3.4 Results</h2>
                        <p>FrogWild is significantly faster and uses less network
                            and CPU compared to GraphLab PR. Let us start with
                            the Twitter graph and the AWS system. In Figure 1(a) we
                            see that, while GraphLab PR takes about 7.5 seconds per
                            iteration (for 12 nodes), FrogWild takes less than 1 sec,
                            achieving more than a 7x speedup. Reducing the value of ps
                            decreases the running time. We see a similar picture when
                            we study the total running time of the algorithms in Figure
                            1(b)).
                            We plot network performance in Figure 1(c). We get a
                            1000x improvement compared to the exact GraphLab PR,
                            and more than 10x with respect to doing 1 or 2 iterations
                            of GraphLab PR. In Figure 1(d) we can see that the total
                            CPU usage reported by the engine is also much lower for
                            FrogWild.
                            We now turn to compare the approximation metrics for
                            the PageRank algorithm. For various k, we check the two
                            accuracy metrics: Mass captured (Figure 2(a)) and the Ex-
                            act identification (Figure 2(b)). Mass captured – is the total</p>
                    </td>
                </tr>
            </tbody></table>
        </div>
          
    
    
</body></html>