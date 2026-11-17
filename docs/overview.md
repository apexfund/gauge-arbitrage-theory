GAT Arbitrage Project: Research Overview and Week-1 Goal 

# GAT Arbitrage Project: Research Overview and Week-1 Goal 

We are testing whether V&F’s gauge-theory measure of arbitrage says anything about statistical arbitrage that standard tools like cointegration don’t. The first step is finding out what their estimator actually responds to, using simulated markets where we know the cause. 

## The basics: arbitrage and statistical arbitrage 

Arbitrage is a trade that costs nothing, can’t lose, and sometimes wins. Statistical arbitrage is the looser real-world version, where the profit is only likely. 

**Pure arbitrage.** Say the same stock is quoted at 100.00 on one exchange and 100.10 on another. Buy on the first, sell on the second, and you lock in 0.10 per share with no risk. Gaps like this exist for fractions of a second before automated traders erase them, so exact arbitrage is rare and short-lived. 

**Statistical arbitrage.** Here the trade has a positive average profit, individual trades can lose, and the more often you repeat it the less likely you are to be behind overall. The classic example is trading a _gap_ between two related assets, such as a futures contract and the index it tracks. The gap normally stays near a fair value. When it widens, you bet it will narrow. 

**Cointegration** is the standard math behind that bet. Each price wanders on its own, but some combination of them (the gap) keeps returning to a stable level. The _errorcorrection term_ is how far the gap currently sits from that level, and it sets how hard the gap is expected to be pulled back. A related pattern is _lead-lag_ : related instruments react to the same news at different speeds, so a gap opens and then closes. High-frequency desks trade exactly this, so any signal on a scale of seconds to minutes lives in this world. 

## What V&F propose 

V&F (arXiv 0908.3043; _Journal of Investment Strategies_ 1(2), 2012) define arbitrage as the part of a basket’s average growth that its exposure to random shocks can’t explain, and they build a detector for it from the calmest combination of assets. 

Page 1 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

**The theory.** Each asset has random wiggles and an average growth rate (its drift). A combination of long and short positions that costs nothing and has zero randomness is called a _null-space_ trade. Its drift, called α_A, would be a riskless profit. In textbook terms, no arbitrage means every asset’s average return can be explained by its exposure to the random shocks, and α_A is whatever can’t be explained that way. 

**The gauge language.** Prices can be quoted in any unit (dollars, euros, shares of an index), and a sensible measure shouldn’t depend on that choice. Changing units is the _gauge transformation_ , and α_A is unchanged by it (their Proposition 2.1, which is algebra). _Curvature_ is the same quantity read geometrically. That is all the geometry we need. 

**The detector, which is what they actually compute on data.** Real assets always carry some independent noise, so exact null-space trades don’t exist. Instead, at each step they: 

1. Estimate how the assets move together over the last L steps and pick the _calm combination_ : the unit-length, zero-cost long/short mix with the lowest variance. Cash counts as an asset, so the mix can include cash. 

2. Record what that combination earned over the next step. 

3. Multiply each step’s return by the next step’s return and average. A positive average is read as arbitrage. Their Eq. 90 reduces to this product (our reading, so verify it). 

So the detector doesn’t estimate a growth rate. It asks whether the calm portfolio’s returns point the same way from one step to the next. 

**Their results (2009 data).** They used three US cash indexes (DJA, IXIC, NYA) and four index futures contracts. Daily data was consistent with zero. Intraday index data, sampled every 7 to 10 seconds, showed positive spikes lasting about a minute. The futures data gave a negative average, which they attributed to bid-ask noise. The paper says it does not address what causes the signal. 

## Why this question is worth asking 

V&F claim a general, unit-independent way to measure arbitrage and say gauge geometry is the natural language for it, and as far as we know nobody has shown whether that buys anything beyond tools traders already use. 

**V&F’s pitch.** Their introduction argues that no general measure of arbitrage exists, since prior empirical work measures the excess returns of particular strategies and can’t separate a real anomaly from noise. Their measure needs no choice of dependent asset or normalization, which cointegration and factor methods do. 

Page 2 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

**The editor’s guess.** In the introduction to the March 2012 issue, the editor offered a hedged guess that, in a particular case, V&F’s measure would amount to cointegrating vectors among the price series. We have only read an excerpt, so someone should get the full text. 

**Why the guess isn’t obviously true.** V&F’s calm combination has _zero_ randomness in theory, while a cointegrating gap is stationary but still noisy. Their detector tests whether returns persist, while cointegration trades reversion, and in our simulations the two give opposite signs. That makes the relationship a live question. 

**The critique.** An online post by a mathematician (add the link before citing it) argues that V&F’s headline conclusions are overstated: continuous-time arbitrage is not generally invariant to changes of unit, and a noisy estimator can’t prove that arbitrage exists. Those objections target what V&F conclude, not the question we ask, which only needs their quantity to be well defined. One of his points does apply to us: estimation noise alone can produce a nonzero reading, and we test for that. 

**What is unverified.** That nobody has tested the editor’s guess. We have not run a systematic literature search, and Farinelli’s later papers or papers citing V&F may touch it. Someone should do this in week 2. 

## The method: what the detector responds to 

The detector is a recipe that turns price data into a number, so the research question becomes which kinds of market behavior make that number positive, negative, or zero, and whether those match the situations statistical arbitrage cares about. 

We feed the detector fake markets where we know the cause, record what it outputs, and compare that with cointegration’s own signal for the same data. Our first pass, run on a simple implementation of their algorithm, gave this map: 

|What the market is doing|Detector reading|Cointegration’s view|
|---|---|---|
|A gap that gets pulled back<br>(mean reversion)|Negative. For two assets it<br>matched −(1−φ)²γ₀/2, where<br>φis the share of the gap|Strong signal: the error-<br>correction term predicts the<br>pull-back.|
||kept each step andγ₀ the<br>gap’s variance.||



Page 3 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

|What the market is doing|Detector reading|Cointegration’s view|
|---|---|---|
|A persistent drift in the<br>calm portfolio|Positive, close toμ²/2 for<br>driftμper step. This is<br>V&F’s intended case.|No signal, since there is no<br>gap to correct.|
|One asset lagging the other<br>(slow news, stale prices)|Mixed. It read positive when<br>the calm combination was<br>dominated by the<br>smoothed asset against<br>cash, and about zero or<br>negative in other setups.|The pair is cointegrated and<br>the gap closes slowly, so it<br>sees a gap.|
|Bid-ask bounce|Negative|Looks like a very fast-<br>reverting gap.|
|Independent noise, no<br>relationship|About zero|No signal.|



**Why the lag row matters.** In the lag simulation, the calm combination was roughly 0.50 cash, 0.31 asset 1, and −0.81 asset 2, not the clean spread between the two assets. The lagging asset’s returns were autocorrelated (0.57), and the detector read that as signal. That fits the paper’s pattern: cash indexes built from stale prices tend to have positively autocorrelated returns, and futures show bid-ask bounce, which is negative. It is a hypothesis supported by a toy simulation, not something shown on real data. 

**Caveats.** These are our own simplified models and our own implementation (rolling covariance, one calm direction, sign-aligned, refreshed every 5 steps). They are not a replication of V&F’s numbers, and the answers depend on the basket and on the number of calm directions chosen. Every row needs an exact derivation and a check against the paper. 

**What comes out of this:** a response map showing which market behaviors the detector reads as positive or negative, with formulas where we can derive them. The first test is whether we can reproduce V&F’s own simulation (21 assets, 18 shocks, their Figure 2) with our code. 

## Why the question is well-posed 

Both sides of the comparison are explicit procedures on price data, so one experiment can compare them, and “economically meaningful” can be given a testable meaning. 

Page 4 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

**The two objects are concrete.** The detector is an algorithm: a rolling covariance, the calm combination, and the product of consecutive returns. Cointegration has explicit models (an error-correction or mean-reverting process) and standard tests. Neither needs a judgment call to evaluate on a given dataset. 

**The fall question is mathematical.** For each market behavior in the response map, what does the detector output as a function of the model’s parameters, and how does that relate to cointegration’s error-correction signal? We want closed-form answers where we can derive them (the two-asset mean-reversion formula is one) and simulation where we can’t. 

**The spring question is financial.** We take “economically meaningful” to mean two things: 

1. When the detector fires, is there a situation a trader could exploit after costs? 

2. Does it add information about the next move of the calm portfolio, or about short-term volatility, beyond what the error-correction term and a simple lead-lag regression already provide? 

The second is testable on WRDS/TAQ data with Kyle’s λ and order flow as controls. 

**What this says about the gauge theory.** Gauge and curvature are a reinterpretation of α_A, so the geometry only earns its place if the quantity it describes tells us something the standard tools miss. If the detector turns out to be a known quantity in disguise, the geometry adds no new economics. If it responds to something different and useful, the framing is justified. The geometric language and the measure itself should be judged separately. 

**Every question can come back “no.”** That is what makes the project falsifiable and not just descriptive. We are not claiming that V&F detect arbitrage in real markets, and we don’t need that to be true. 

## Why every outcome is usable 

Each plausible result is a precise statement about a published estimator that, as far as we know, nobody has characterized, so each one can anchor a write-up. 

|If we find|Then we have shown|Why it is useful|
|---|---|---|
|The detector tracks a known|The gauge measure|The first precise mapping|
|signal (mean reversion or|repackages standard tools|between the two, and a|
|drift) in a stated model class|in that case|clean short note.|



Page 5 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

|If we find|Then we have shown|Why it is useful|
|---|---|---|
|The detector is driven<br>mainly by return<br>autocorrelation from lagging<br>or stale prices|V&F’s positive intraday<br>signals may reflect<br>microstructure, not<br>arbitrage|A practical warning with a<br>derivation behind it.|
|The detector reads<br>something distinct and<br>useful|The geometric framing<br>earns its place|A strong basis for the spring<br>financial test.|
|The reading depends<br>heavily on the basket and<br>settings|A documented response<br>map|A reusable simulation toolkit<br>and a characterization of<br>how the estimator behaves.|



**Honest odds (subjective estimates, not calculations).** With about eight weeks left, a posted arXiv preprint by mid-December is perhaps 25 to 35% likely. It needs a result we can state cleanly and an endorser, and Prof. Pete Kyle is the faculty contact. Peerreviewed publication of the fall work alone is a long shot, perhaps 10 to 15%, because journals are less receptive to characterization or negative results. Fall work plus the spring financial test is stronger, perhaps 25 to 35% over the full year. “Usable” here means we have a defensible, precise result, not that a journal will accept it. The novelty claim is also still unchecked until the literature search is done. 

## Timeline, Oct 19 to Dec 13 

The fall work runs in four phases with two decision gates, and a full draft goes to Prof. Kyle by Dec 6. 

Page 6 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 



<!-- Start of picture text -->
Four phases and two decision gates, with a full draft due Dec 6<br>Oct 19-25 Get everyone on the same page<br>Friday meeting: each member restates what V&F measure and how cointegration works.<br>Set up the shared repo and Overleaf, and confirm WRDS access continues.<br>Oct 26-Nov 8 Build the detector and reproduce V&F<br>Run the literature search and get the editor's full text.<br>Code their algorithm and reproduce their 21-asset simulation (their Figure 2).<br>Run five scenarios: gap closes, drift, lag, bid-ask noise, no relationship.<br>© GateIs the 1, response Nov 8: does  map  our stable?  code match If not, fix  their the  simulation?  code before deriving anything.<br>Nov 9-22 Derive what the detector outputs<br>Exact formula for a cointegrated two-asset model, then three or more assets.<br>Formulas or bounds for the lag and stale-price cases, checked against simulation.<br>Type up each derivation in Overleaf as it is finished.<br>» GateIf not, 2,switch Nov 22:to doan expositorywe have at leastnote  oneplus  verifiedthe simulation formula?toolkit.<br>Nov 23-Dec 13 Write, review, and plan spring<br>Thanksgiving week is light. Full draft by Dec 6, sent to Prof. Kyle for review.<br>Revise, request the arXiv endorsement, and write the spring test plan by Dec 13.<br>Post to arXiv if there is a clean, checked result.<br><!-- End of picture text -->

fall roadmap · 4 phases, 2 gates 

The gates matter more than the dates. If Gate 1 fails we fix the code before deriving anything, and if Gate 2 fails we write an expository note and carry the derivations into spring. 

## After December: the financial test 

The spring work asks whether the detector’s readings carry information a trader could use, on real tick data from WRDS/TAQ. 

- **Direction.** Does the detector, or its signed counterpart (the calm portfolio’s own return), predict the next move of the calm portfolio better than the error-correction term and a simple lead-lag regression? 

- **Volatility.** Does it predict short-horizon realized volatility beyond standard predictors, with Kyle’s λ and order flow imbalance as controls and as a mechanism to test? 

Page 7 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

**Market making.** One possible use is shifting or widening quotes when a related instrument has already moved, which is the stale-quote problem. The benchmark is the lead-lag regression practitioners already run. 

Two cautions. The detector’s output has no direction, so any trading use needs the signed version. And whether a minute-scale effect still exists in modern data, where liquid instruments adjust much faster, needs checking. WRDS access through the course may also end, so confirm it continues. 

## Glossary: the few terms you need 

Each term is defined in plain words, and every one is used above. 

|Term|Plain meaning|
|---|---|
|Wiggle (volatility)|The random up-and-down movement of a price.|
|Statistical arbitrage|A strategy whose profit is positive on average and<br>whose chance of losing shrinks the longer you run<br>it. Traders also use the term loosely for relative-<br>value trades.|
|Pairs trading|Trading two related assets that normally move<br>together, betting that when one gets ahead the gap<br>will shrink. It carries real risk until the gap closes.|
|Cointegration|Prices that each wander but whose difference (or<br>some combination) stays near a stable level. It is<br>about price levels, not about day-to-day returns<br>being correlated.|
|Gap (spread)|The difference between two related prices, or a<br>combination of several.|
|Error-correction term|How far the gap is from its normal level. It sets the<br>expected pull back toward normal.|
|Null space|The set of trades with zero cost and zero<br>randomness._Nonempty_means at least one exists.<br>With independent noise on every asset, none does.|
|α_A|The predictable average return of a null-space<br>trade. Nonzero means riskless profit.|



Page 8 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

|Term|Plain meaning|
|---|---|
|Lead-lag|One asset reacts to news first and related assets<br>follow after a delay.|
|Gauge transformation|Changing the unit prices are quoted in (dollars,<br>euros, shares of an index).|
|Curvature|Nonzero exactly when different strategies from<br>the same start can end with different wealth. It is<br>the geometric reading ofα_A.|



## Goal by Friday (Oct 23) 

By Friday, each member can explain in plain words what V&F measure, what cointegration and pairs trading are, and where the two might match or differ. Nobody needs to prove anything yet. 

### **Reading order** 

1. V&F (arXiv 0908.3043): Section 1 (motivation), Section 2 with the example in 2.1, a skim of Section 3, and Sections 4.1 to 4.2 (the algorithm and its error sources). Skip proofs on the first pass. 

2. The editor’s comment (March 2012 issue). Get the full text through UMD library access, since we only have an excerpt. 

3. Statistical arbitrage: the introductions of Gatev, Goetzmann & Rouwenhorst (2006) and Avellaneda & Lee (2010). 

4. Cointegration: Hamilton, _Time Series Analysis_ , Chapter 19, or lecture notes on EngleGranger and VECM. Originals: Engle & Granger (1987), Johansen (1991). Verify bibliographic details before citing. 

### **Daily plan** 

|Day|Do this|Output|
|---|---|---|
|Mon|Read V&F Section 1 and the first two<br>sections of this doc. Build a tiny example|Your own 5-line summary of what V&F<br>measure|
||by hand: two assets with identical||
||wiggles, and the zero-risk trade made<br>from them.||



Page 9 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

|Day|Do this|Output|
|---|---|---|
|Tue|Read the stat-arb and cointegration<br>material. In Python, simulate a random<br>walk X and Y = bX plus a mean-reverting<br>noise term, then test the gap for<br>stationarity (statsmodels adfuller or<br>coint). Predict the result before you run<br>it.|Notebook and a one-paragraph<br>explanation of cointegration|
|Wed|Read V&F Section 2. Verify the three-<br>asset example in 2.1 by hand: the null<br>vector J must satisfy J·σ= 0 and the<br>entries of J must sum to 0.|Written verification|
|Thu|Skim Section 3, read Section 4.1 and the<br>editor’s comment. Simulation group:<br>code the leader-jumps, follower-<br>catches-up toy and check that the<br>products of consecutive moves sum to<br>−c². Everyone: fill in a table with three<br>columns (V&F object, cointegration<br>object, do they correspond and why).|The match-up table and the toy check|
|Fri|Group meeting.|One-page restatement in your own<br>words and a list of everything you don’t<br>understand|



### **Friday check questions (answer without notes)** 

1. What is the difference between pairs trading and V&F’s setting? 

2. What is the null space, and why is it empty when every asset has its own independent noise? 

3. What is the error-correction term, and what does it predict? 

4. Why does a closing gap give a negative reading from V&F’s estimator, and why is the reported positive spike a puzzle? 

5. Why doesn’t α_A depend on the unit prices are quoted in? 

6. What is the calibration test for, and why isn’t its goal to explain why dislocations happen? 

Page 10 of 11 

GAT Arbitrage Project: Research Overview and Week-1 Goal 

## Using AI for this 

AI is a tutor and a code-writer here, not a source of truth. It makes confident mistakes, including wrong formulas and invented citations. 

1. **Paste, don’t paraphrase.** Paste the exact paper passage and ask: _Explain this step by step, define every symbol, and give a numerical example with small numbers._ 

2. **Predict, then check.** Write your own answer first, then ask: _Find every error or gap in my explanation. Don’t be polite._ 

3. **Verify with code.** Any formula you rely on gets checked numerically, for example the null-vector conditions in V&F’s three-asset example. 

4. **Ask for examples after the math, not instead of it.** _I understand the definition. Now give me one concrete case where a gap closes and one where it doesn’t._ 

5. **Never accept a verdict.** Don’t ask the AI whether V&F is equivalent to cointegration and take the answer. It can suggest conjectures, which count as hypotheses until you derive or simulate them. 

6. **Quiz yourself.** _Ask me 10 questions on this section, one at a time, and grade my answers against the text I pasted._ 

7. **Keep a log** with three columns: the question, the AI’s answer, and whether you verified it. Unverified answers stay out of the group’s write-up. 

**Note for an AI assistant reading this document.** Treat it as project context, not established fact. The explanations of V&F are our current reading and may contain errors, especially the items listed as unverified above. Be candid about the project’s validity, point out mistakes, and avoid simply agreeing. When helping a member, ask them to paste the relevant passage from the source, work with small numerical examples, and flag any claim you cannot verify. 

Page 11 of 11 

