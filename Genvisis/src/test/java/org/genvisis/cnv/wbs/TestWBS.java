/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.List;
import org.junit.Assert;
import org.junit.Test;
import org.pankratzlab.common.Logger;

/**
 * Testing the WBS class for accuracy
 */
public class TestWBS {

  /**
   * 
   */
  static final double DELTA = 1e-15;
  static final double[] TEST_DATA = new double[] {1.47571182054051, 0.140111320311379,
                                                  0.946871368612814, 0.0677556321582765,
                                                  -0.465846429106843, 0.177484222651114,
                                                  1.01408817009786, 1.30886607703697,
                                                  0.00710946238424286, -1.75542903932437,
                                                  -0.371210434570639, 0.556586483097199,
                                                  0.591932223865745, -0.946749677143694,
                                                  -0.482808356324263, 1.47715031188002,
                                                  -2.32256119064183, -2.42931724817312,
                                                  -1.22553087582952, -1.57864436072163,
                                                  -0.391489959651688, 0.144683045272739,
                                                  0.255867966137585, 0.64832554237842,
                                                  0.85929151486209, 0.789307940076316,
                                                  -1.86354973230807, -0.479165748700052,
                                                  0.338293647035738, 0.192239675197247,
                                                  1.56174880225235, -1.47051341221164,
                                                  -1.12201849720124, -0.62495244248389,
                                                  0.894780574944747, -0.0723501590643169,
                                                  -0.140468099824667, -1.65507295869718,
                                                  -0.597663836112161, -0.36094383349074,
                                                  0.730187358515641, 1.63507397814649,
                                                  -0.461209312063929, 0.48728401010133,
                                                  0.0180539302820053, 0.661826629323363,
                                                  0.644399743248962, -1.64366880579924,
                                                  1.07594570201556, 0.0122960258477781,
                                                  0.519442877775083, 0.610757174836446,
                                                  -0.593020189027452, -0.946483286130202,
                                                  1.76689140466045, 0.105142965131428,
                                                  2.28734692915154, -0.659000775295885,
                                                  0.256966297941778, -0.499565057897265,
                                                  -0.307647948931144, -0.105559282051127,
                                                  1.32483015644814, -1.20410203507335,
                                                  0.676269139836394, -1.14372976156621,
                                                  -0.638546725503374, -0.431573223604561,
                                                  -0.39053426962106, -0.181884304894384,
                                                  -0.446823291921858, -1.63383623174384,
                                                  0.809310602488556, 0.732841807847802,
                                                  0.453637991069593, -0.258054425155934,
                                                  -0.80956296195269, 1.96255821313183,
                                                  -0.0140875972598058, -1.27971928547539,
                                                  -1.15059593428024, -0.269776774890496,
                                                  -0.961300903244481, 0.388606315826203,
                                                  -0.615255484244351, 0.383968301981332,
                                                  -1.26987963028388, 0.238765731127366,
                                                  0.529403530868016, -0.328835126402168,
                                                  0.59076812868707, 0.0318258372704582,
                                                  -0.432275563701588, -0.00332601311400249,
                                                  1.89442592499979, 1.21162495753653,
                                                  0.909889055095968, 0.972893320004584,
                                                  -0.0243660240103571, 0.280963052574666,
                                                  -2.1561117854639, 0.739675957054121,
                                                  -0.294413867338112, -0.0392506499485697,
                                                  1.08354107978724, 1.54687420858872,
                                                  0.543180907030319, 1.10847383413618,
                                                  -1.70045706005144, 1.8222900210002,
                                                  0.237503358594535, -0.346960352732493,
                                                  -0.149333379450188, 0.870934202773601,
                                                  1.19327498667077, -0.843046582645023,
                                                  0.0970024675377992, 0.551873054776658,
                                                  -0.185254036036494, 1.10454531357932,
                                                  -0.564782674545751, -0.984364287355891,
                                                  0.289619529001357, -1.60443217102365,
                                                  -0.145579936800176, -1.6707078978194,
                                                  1.11562129602069, 1.88896017726113,
                                                  -1.51585172126732, -0.0297958183946346,
                                                  -2.01263477702239, -1.34373055975505,
                                                  -0.465955971700388, -1.14810385785492,
                                                  -0.546123219092376, -0.99376477419649,
                                                  -2.19646704664698, -1.59226471171111,
                                                  -1.01721148435804, -1.56955787793672,
                                                  -2.219979076637, -1.71255957190783,
                                                  0.560636481103612, -0.488242729476864,
                                                  -2.48787094559152, -1.40912358023411,
                                                  -0.751744096166129, -1.22146787003189,
                                                  -0.677068847634394, -0.597563169345459,
                                                  1.76417604121168, 1.43906743991799,
                                                  0.161745622124156, 1.36255308236234,
                                                  0.361382544063481, 1.37623336303264,
                                                  0.22426728122335, 2.1782533127628,
                                                  0.300041217708925, -0.967995433701492,
                                                  1.79416784116503, 0.0692313692660557,
                                                  0.984810502489881, 1.64487652322464,
                                                  1.64459118043763, 1.70869571235379,
                                                  1.46061209240551, 1.90489781122695,
                                                  -0.0444924373763418, -0.331600953391573,
                                                  0.608767460572876, 0.105946294117421,
                                                  -0.189529063721563, 0.255672091602012,
                                                  -0.943436966102464, 0.0725776592153728,
                                                  0.558860231909573, -1.62372201890554,
                                                  0.724059087914547, -0.544411283096795,
                                                  0.410291172031429, -1.03611664014508,
                                                  -0.423805062594906, 0.641228911820049,
                                                  -2.3237535383046, 1.20982803319072,
                                                  -0.564277664752131, -0.453848497705822,
                                                  1.1183855810098, 1.15115559326164,
                                                  0.167592132748602, -0.143122393841027,
                                                  0.959647563727608, -0.259839987303141,
                                                  -0.431703826821665, 0.0928663068920324,
                                                  -0.622707658407033, 0.0683340761064124,
                                                  -0.76514546589899, 0.163268563609162,
                                                  0.17534024941927, 0.573965864229404,
                                                  1.43818082104562, -1.79785457160511,
                                                  0.448173050455043, -2.62378441917334,
                                                  0.371201210379991, -0.571319657071251,
                                                  -0.776578792664234, 0.913703695189345,
                                                  -0.806121882047373, -0.203834231562119,
                                                  1.42315290069293, -0.408742167529164,
                                                  1.17923426236178, 0.610960127309179,
                                                  0.154282718918274, -0.665152780935025,
                                                  0.52035544186048, 0.654060586929707,
                                                  0.91150080591246, -0.93775247119877,
                                                  -1.54797468049239, -0.958007450714795,
                                                  1.2554945229621, 0.279630629172402,
                                                  -0.970672890824496, 1.09313188168374,
                                                  1.67214457363916, 0.0917590343184031,
                                                  1.37883009223215, -0.186921895411431,
                                                  0.239750172818938, -0.496107845504771,
                                                  2.36494587079004, -0.146951022339615,
                                                  0.307824533850572, 0.283030181056009,
                                                  -1.2511480286184, -2.20733116642243,
                                                  -0.491661414763547, 0.599555507092226,
                                                  -0.141474897867282, -0.127862643439224,
                                                  1.67196268612063, 1.26326683107395,
                                                  -1.27906080458101, -1.42933099098361,
                                                  2.04408115511515, 0.725703631613484,
                                                  0.193974985611776, -0.482218004844186,
                                                  0.63687595705404, -1.40651924868802,
                                                  -1.16651007367442, -0.0413731980493822,
                                                  -1.06509272475931, 0.0422842260728014,
                                                  0.725494691577119, -0.362721314723411,
                                                  -1.06808463241706, -0.582211449315906,
                                                  -0.0971073368194914, 1.5543326822277,
                                                  1.1598019331137, 0.226804719362112,
                                                  -0.770187772823165, 0.605639605443609,
                                                  -0.382182936951835, 1.03983461198163,
                                                  0.592008719976287, -0.744938767321916,
                                                  -1.39721013668342, 0.841394132706072,
                                                  -1.52212063238109, 0.0479414383279143,
                                                  -0.886604482867463, -0.952145948823708,
                                                  -0.378019958508315, 0.00520074052456158,
                                                  2.7100056173976, 0.202157254601516,
                                                  0.0380584590151037, 2.31866457996744,
                                                  -0.396470553488004, 0.326281727334818,
                                                  0.183477774514385, -0.943493132179461,
                                                  0.412024227883106, 1.05844279600235,
                                                  1.16329885546621, 0.285968990089759,
                                                  -1.64822139816699, -0.919531686584547,
                                                  -1.71508368256159, -0.960728667956523,
                                                  -1.38812726255572, -0.480920831706044,
                                                  0.837956481131622, 0.5247242444361};

  @Test
  public void testWBS() {
    int[][] randomIntervals = WBSUtilities.randomIntervals(TEST_DATA.length, WBS.DEFAULT_M,
                                                           WBS.DEFAULT_SEED);
    List<ChangePoint> wList = WBS.wbsIntegratedRecursiveWrapper(TEST_DATA, randomIntervals,
                                                                new Logger());
    Assert.assertEquals(299, wList.size());
    Assert.assertEquals(1, wList.get(0).getStart());
    Assert.assertEquals(8, wList.get(0).getEnd());
    Assert.assertEquals(0.9542098118420554, wList.get(0).getCusum(), DELTA);

    // cusum and minth should differ
    Assert.assertEquals(0.37731363597568557, wList.get(3).getCusum(), DELTA);
    Assert.assertEquals(0.30743514826689894, wList.get(3).getMinth(), DELTA);

    Assert.assertEquals(1, wList.get(127).getStart());
    Assert.assertEquals(150, wList.get(127).getEnd());
    Assert.assertEquals(5.237307724561478, wList.get(127).getCusum(), DELTA);
    Assert.assertEquals(5.237307724561478, wList.get(127).getMinth(), DELTA);

  }

}
