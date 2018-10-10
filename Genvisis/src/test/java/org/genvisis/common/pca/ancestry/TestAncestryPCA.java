/**
 * 
 */
package org.genvisis.common.pca.ancestry;

import static org.junit.Assert.assertEquals;
import java.util.HashMap;
import java.util.Map;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.genvisis.cnv.ejml.matrix.MatrixDataLoading;
import org.genvisis.cnv.ejml.matrix.NamedRealMatrix;
import org.genvisis.cnv.gwas.pca.ancestry.AncestryPCA;
import org.junit.Test;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;

/**
 * A few basic tests for {@link AncestryPCA}.<br>
 * Note, the initial (non-extrapolated) PC output we test against was computed using the following
 * "AncestryPCA" implementation in R
 * 
 * <pre>
 * input = data.frame(
  M1 = c(0, 1, 2, 2, 2, 1, 2, 1, 0, 0),
  M2 = c(0, 1, 2, 1, 2, 0, 2, 1, 0, 0),
  M3 = c(1, 0, 2, 1, 2, 1, 2, 0, 0, 0),
  M4 = c(0, 1, 2, 2, 0, 1, 0, 1, 1, 1),
  M5 = c(0, 0, 0, 0, 1, 2, 1, 2, 2, 2),
  M6 = c(0, 1, 2, 2, 2, 1, 2, 1, 0, 1),
  M7 = c(0, 1, 0, 1, 2, 0, 2, 1, 2, 0),
  M8 = c(1, 0, 2, 0, 2, 1, 2, 1, 0, 0),
  M9 = c(0, 1, 2, 2, 2, 1, 2, 1, 1, 1),
  M10 = c(0, 0, 0, 0, 1, 1, 1, 2, 2, 2),
  M11 = c(0, 1, 0, 0, 1, 1, 1, 2, 2, NA)
)


for (marker in colnames(input)) {
  tmp = input[, marker]
  mean = mean(tmp, na.rm = TRUE)
  sum = sum(tmp, na.rm = TRUE)
  count = length(which(!is.na(tmp)))
  
  possible = 2 + 2 * count
  pi = (1 + sum) / possible
  norm = sqrt(pi * (1 - pi))
  input[, marker] = input[, marker] - mean
  input[, marker] = input[, marker] / norm
  input[, marker][is.na(input[, marker])] <- 0
  print(paste0(" marker=", marker, " mean=", mean, " sum=", sum, " count=", count))
  
}

input = as.data.frame(t(input))
colnames(input) = paste0("SAMPLE_", seq(0:9))
rownames(input) = paste0("MARKER_", seq(0:10))

# PCs
svd = svd(x = input)
v = as.data.frame(svd$v)
v = v[, c(1:5)]
colnames(v) = paste0("PC", seq(1:5))
print(v, digits = 16)

# Loadings
u = as.data.frame(svd$u)
u = u[, c(1:5)]
colnames(u) = paste0("PC", seq(1:5))
print(u, digits = 16)
 * </pre>
 * 
 * With the following sessionInfo()
 * 
 * <pre>
 * R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin17.5.0 (64-bit)
Running under: macOS High Sierra 10.13.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_3.5.0 tools_3.5.0    yaml_2.1.19
 * </pre>
 */
public class TestAncestryPCA {

  private static final double[][] INPUT = new double[][] {{0, 1, 2, 2, 2, 1, 2, 1, 0, 0},
                                                          {0, 1, 2, 1, 2, 0, 2, 1, 0, 0},
                                                          {1, 0, 2, 1, 2, 1, 2, 0, 0, 0},
                                                          {0, 1, 2, 2, 0, 1, 0, 1, 1, 1},
                                                          {0, 0, 0, 0, 1, 2, 1, 2, 2, 2},
                                                          {0, 1, 2, 2, 2, 1, 2, 1, 0, 1},
                                                          {0, 1, 0, 1, 2, 0, 2, 1, 2, 0},
                                                          {1, 0, 2, 0, 2, 1, 2, 1, 0, 0},
                                                          {0, 1, 2, 2, 2, 1, 2, 1, 1, 1},
                                                          {0, 0, 0, 0, 1, 1, 1, 2, 2, 2},
                                                          {0, 1, 0, 0, 1, 1, 1, 2, 2, Double.NaN}};

  /**
   * Additional samples for extrapolation
   */
  private static final double[][] ADDITIONAL_SAMPLES = Matrix.transpose(new double[][] {{0, 1, 2, 2,
                                                                                         1, 2, 0, 0,
                                                                                         Double.NaN,
                                                                                         2, 1},
                                                                                        {0,
                                                                                         Double.NaN,
                                                                                         2, 2, 1, 2,
                                                                                         0, 0, 1, 0,
                                                                                         0}});

  /**
   * Additional markers for extrapolation (includes data for
   * {@link TestAncestryPCA#ADDITIONAL_SAMPLES}
   */
  private static final double[][] ADDITIONAL_MARKERS = new double[][] {{1, 0, 2, 2, 1, 0, 1, 0,
                                                                        Double.NaN, 2, 1, 1},
                                                                       {0, 0, 2, 2, 1, 2, 0, 0, 1,
                                                                        2, 1, 1}};

  /**
   * The output PCs generated by the R script above, note that the +/- sign of an entire PC can be
   * flipped without consequence, but must be consistent
   */
  private static final double[][] OUTPUT_PCS = new double[][] {{0.16119448641633580,
                                                                -0.52666965438145663,
                                                                0.68436660140910732,
                                                                -0.1400052060948393,
                                                                -0.02562962410622701},
                                                               {0.04710839548459887,
                                                                -0.26975880222591397,
                                                                -0.18400549682603928,
                                                                -0.5143397982904894,
                                                                0.37335365918393359},
                                                               {-0.46994363022424829,
                                                                -0.25035636492945640,
                                                                -0.11407490849750031,
                                                                0.4026567795625129,
                                                                0.37514833479108606},
                                                               {-0.24120269509174444,
                                                                -0.29829188974362808,
                                                                -0.56421111413995795,
                                                                -0.1782070609993075,
                                                                -0.42794580875879068},
                                                               {-0.35959732246293818,
                                                                0.39209515775686077,
                                                                0.23651066829344944,
                                                                -0.1235838682165552,
                                                                -0.15927078858573421},
                                                               {0.15470686408036172,
                                                                -0.03236721991926313,
                                                                0.08180557677897178,
                                                                0.4837095576706957,
                                                                -0.19399637599744360},
                                                               {-0.35959732246293818,
                                                                0.39209515775686077,
                                                                0.23651066829344941,
                                                                -0.1235838682165554,
                                                                -0.15927078858573404},
                                                               {0.21923139703354966,
                                                                0.32901202310356492,
                                                                -0.12114062376593708,
                                                                0.1379004388734949,
                                                                0.60067258009343916},
                                                               {0.45952972872598624,
                                                                0.28461230336264887,
                                                                -0.11225489288749005,
                                                                -0.3164303875712338,
                                                                -0.10927540648330344},
                                                               {0.38857009850103708,
                                                                -0.02037071078021766,
                                                                -0.14350647865805335,
                                                                0.3718834132822771,
                                                                -0.27378578155122557}};

  /**
   * The output (marker) loadings generated by the R script above, same deal with the sign
   */
  private static final double[][] OUTPUT_LOADINGS = new double[][] {{-0.41648070178497237,
                                                                     0.12071744077230184,
                                                                     -0.22284269571024129,
                                                                     0.02519086503736782,
                                                                     0.02806017855167750},
                                                                    {-0.40171462398676744,
                                                                     0.20087550582841796,
                                                                     -0.05324867070269729,
                                                                     -0.09926901036725208,
                                                                     0.49528871838788274},
                                                                    {-0.39325375584602212,
                                                                     0.05099021241131581,
                                                                     0.32338295618177038,
                                                                     0.19411319606361638,
                                                                     -0.40143191987660615},
                                                                    {-0.02614424518637623,
                                                                     -0.19543069738869459,
                                                                     -0.64534776546813855,
                                                                     0.24917138492628391,
                                                                     0.21889192354919262},
                                                                    {0.29324268465020825,
                                                                     0.46012852857952741,
                                                                     -0.04102212867514955,
                                                                     0.44910077042381880,
                                                                     -0.20297640700940972},
                                                                    {-0.35460685914750512,
                                                                     0.11725414618613965,
                                                                     -0.27677528350807412,
                                                                     0.17894489446890002,
                                                                     -0.17988171012405496},
                                                                    {-0.08436578140932870,
                                                                     0.46024869146156622,
                                                                     -0.05196898484101864,
                                                                     -0.68517256403146065,
                                                                     -0.23254934387938561},
                                                                    {-0.31465105851632952,
                                                                     0.20306062688507953,
                                                                     0.47914839004404269,
                                                                     0.32289348583962268,
                                                                     0.37130850463803683},
                                                                    {-0.28123405370432347,
                                                                     0.19125991764860589,
                                                                     -0.32373519996797462,
                                                                     0.04945916552359558,
                                                                     -0.26882656461625437},
                                                                    {0.26805131563644713,
                                                                     0.46988818365725010,
                                                                     -0.06995218919101083,
                                                                     0.25390778362487548,
                                                                     -0.05808233374150238},
                                                                    {0.20256573437165667,
                                                                     0.41033965976331477,
                                                                     -0.07865317001410922,
                                                                     -0.12407649800349710,
                                                                     0.45132373923309654}};

  /**
   * Shooting to be within this tolerance of the R method above
   */
  private static final double DELTA = 1e-14;

  /**
   * Small implementation of {@link MatrixDataLoading} for testing
   */
  private static class InputProvider implements MatrixDataLoading {

    private final double[][] data;

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.matrix.MatrixDataLoading#getData()
     */
    /**
     * @param data
     */
    public InputProvider(double[][] data) {
      super();
      this.data = data;
    }

    @Override
    public NamedRealMatrix getData() {

      DenseMatrix64F m = new DenseMatrix64F(data.length, data[0].length);

      for (int i = 0; i < data.length; i++) {
        for (int j = 0; j < data[0].length; j++) {
          m.add(i, j, data[i][j]);
        }
      }

      Map<String, Integer> mockMarkers = new HashMap<>();
      Map<String, Integer> mockSamples = new HashMap<>();

      for (int i = 0; i < m.getNumCols(); i++) {
        mockSamples.put("Sample_" + i, i);
      }
      for (int i = 0; i < m.getNumRows(); i++) {
        mockMarkers.put("Marker_" + i, i);
      }
      return new NamedRealMatrix(mockMarkers, mockSamples, m);
    }

  }

  @Test
  public void testAncestryPCsAndLoadings() {
    int numberOfPCSamples = INPUT[0].length;
    int numberOfMarkers = INPUT.length;
    int nonPCSamples = ADDITIONAL_SAMPLES[0].length;

    InputProvider inputProvider = new InputProvider(INPUT);
    NamedRealMatrix m = inputProvider.getData();
    assertEquals(10, m.getDenseMatrix().getNumCols());//samples
    assertEquals(numberOfPCSamples, m.getDenseMatrix().getNumCols());
    assertEquals(11, m.getDenseMatrix().getNumRows());//markers
    assertEquals(11, numberOfMarkers);//markers

    Logger log = new Logger();
    AncestryPCA ancestryPCA = AncestryPCA.generatePCs(inputProvider, 5, log);
    //compare the generated PCs to expected
    DenseMatrix64F pcs = ancestryPCA.getSvd().getPCs().getDenseMatrix();
    DenseMatrix64F transposed = new DenseMatrix64F(pcs.numCols, pcs.numRows);
    CommonOps.transpose(pcs, transposed);
    comparePCsToStandard(transposed, OUTPUT_PCS, numberOfPCSamples);
    //compare the generated loadings to expected

    compareLoadingsToStandard(ancestryPCA.getSvd().getLoadings(), OUTPUT_LOADINGS);

    //compare the extrapolated pcs on the same samples to expected
    comparePCsToStandard(AncestryPCA.extrapolatePCs(ancestryPCA, inputProvider, log)
                                    .getDenseMatrix(),
                         OUTPUT_PCS, numberOfPCSamples);

    //    spike in new samples
    double[][] inputWithNewSamples = new double[INPUT.length][numberOfPCSamples + nonPCSamples];
    for (int i = 0; i < inputWithNewSamples.length; i++) {
      inputWithNewSamples[i] = ArrayUtils.concatDubs(INPUT[i], ADDITIONAL_SAMPLES[i]);
    }

    //compare the extrapolated pcs generated with the spiked in samples to expected for the original samples

    compareExtrapolated(numberOfPCSamples, nonPCSamples, 11, log, ancestryPCA, inputWithNewSamples);

    //spike in new markers
    double[][] inputWithNewSamplesNewMarkers = new double[INPUT.length
                                                          + ADDITIONAL_MARKERS.length][numberOfPCSamples
                                                                                       + nonPCSamples];
    for (int i = 0; i < INPUT.length; i++) {
      inputWithNewSamplesNewMarkers[i] = ArrayUtils.concatDubs(INPUT[i], ADDITIONAL_SAMPLES[i]);
    }
    for (int i = 0; i < ADDITIONAL_MARKERS.length; i++) {
      inputWithNewSamplesNewMarkers[i + INPUT.length] = ADDITIONAL_MARKERS[i];
    }

    compareExtrapolated(numberOfPCSamples, nonPCSamples, 13, log, ancestryPCA,
                        inputWithNewSamplesNewMarkers);

  }

  private void compareExtrapolated(int numberOfPCSamples, int nonPCSamples, int numberOfMarkers,
                                   Logger log, AncestryPCA ancestryPCA, double[][] input) {
    InputProvider inputProvider = new InputProvider(input);
    NamedRealMatrix m = inputProvider.getData();
    assertEquals(numberOfPCSamples + nonPCSamples, m.getDenseMatrix().getNumCols());// samples
    assertEquals(numberOfMarkers, m.getDenseMatrix().getNumRows());// number of markers
    DenseMatrix64F pcs = AncestryPCA.extrapolatePCs(ancestryPCA, inputProvider, log)
                                    .getDenseMatrix();
    comparePCsToStandard(pcs, OUTPUT_PCS, numberOfPCSamples);
  }

  private void comparePCsToStandard(DenseMatrix64F pcs, double[][] standardPCs, int numSamples) {
    // ensure we test all PCs we think we have
    assertEquals(5, pcs.numCols);
    for (int component = 0; component < pcs.numCols; component++) {//
      boolean flip = false;

      for (int sample = 0; sample < numSamples; sample++) {

        double standard = standardPCs[sample][component];
        double test = pcs.get(sample, component);
        //since the sign of a PC does not matter, we test if the first row (sample) of the PC should be flipped. 
        //  If it should, we flip all other entries.
        //   We only do this test once per PC since it should be consistent      
        if (sample == 0 && ((standard < 0 && test > 0) || (standard > 0 && test < 0))) {
          flip = true;
        }
        if (flip) {
          test = test * -1;
        }
        assertEquals(standard, test, DELTA);
      }
    }
  }

  void compareLoadingsToStandard(NamedRealMatrix loadings, double[][] standardLoadings) {
    for (int component = 0; component < loadings.getDenseMatrix().numCols; component++) {//
      boolean flip = false;

      for (int marker = 0; marker < loadings.getDenseMatrix().numRows; marker++) {

        double standard = standardLoadings[marker][component];
        double test = loadings.getDenseMatrix().get(marker, component);
        if (marker == 0 && ((standard < 0 && test > 0) || (standard > 0 && test < 0))) {
          flip = true;
        }
        if (flip) {
          test = test * -1;
        }
        assertEquals(standard, test, DELTA);
      }
    }
  }
}
