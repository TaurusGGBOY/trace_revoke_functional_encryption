package algorithm;

import it.unisa.dia.gas.jpbc.Element;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GaussTest {
    @Test
    void twoIdTest() throws IOException {
        int lambda = 256;
        int t = 6;
        int r = 3;
        int id = 2;
        List<Integer> revokeId = new ArrayList<>(Arrays.asList(1, 3));

        TraceRevokeFE fe = new TraceRevokeFE();

        // setup
        TraceRevokeFE.SetupParams setupParams = fe.setup(lambda, t, r);
        TraceRevokeFE.Params params = setupParams.PK.params;
        // init revoke tree
        List<TraceRevokeFE.Theta> gamma = new ArrayList<>();
        for (int revoke : revokeId) {
            gamma.add(params.idTree.gamma.get(revoke));
        }
        TraceRevokeFE.IdTree revokeTree = fe.new IdTree(gamma);

        LocalDateTime startSetup = LocalDateTime.now();
        List<Element> vr = fe.getVR(params, revokeTree);
        LocalDateTime endSetup = LocalDateTime.now();

        for (int revoke : revokeId) {
            Element sum = params.pairing.getZr().newZeroElement();
            TraceRevokeFE.Theta theta = params.idTree.gamma.get(revoke);
            for (int i = 0; i < theta.theta.size(); i++) {
                sum.add(theta.theta.get(i).mul(vr.get(i)));
            }
            assert sum.equals(params.pairing.getZr().newZeroElement());
        }
        Duration setupDuration = Duration.between(startSetup, endSetup);
        System.out.println("getvr: " + setupDuration.toMillis());

    }
}