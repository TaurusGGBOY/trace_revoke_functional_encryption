package algorithm;


import it.unisa.dia.gas.jpbc.Element;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.*;

public class TraceRevokeFETest {
    @Test
    public void testTraceableFE() throws IOException {
        int lambda = 256;
        int t = 6;
        int r = 3;
        int id = 4;
        List<Integer> revokeId = new ArrayList<>(Arrays.asList(1, 2));

        TraceRevokeFE fe = new TraceRevokeFE();

        // setup
        LocalDateTime startSetup = LocalDateTime.now();
        TraceRevokeFE.SetupParams setupParams = fe.setup(lambda, t, r);
        LocalDateTime endSetup = LocalDateTime.now();
        TraceRevokeFE.Params params = setupParams.PK.params;

        // init m
        Element m = params.pairing.getZr().newRandomElement().getImmutable();


        // keyGen
        LocalDateTime startKeyDer = LocalDateTime.now();
        TraceRevokeFE.SKX SKX = fe.keyGen(params, id, setupParams.msk);
        LocalDateTime endKeyDer = LocalDateTime.now();

        // init revoke tree
        List<TraceRevokeFE.Theta> gamma = new ArrayList<>();
        for (int revoke : revokeId) {
            gamma.add(params.idTree.gamma.get(revoke));
        }
        TraceRevokeFE.IdTree revokeTree = fe.new IdTree(gamma);

        // encrypt
        LocalDateTime startEncrypt = LocalDateTime.now();
        TraceRevokeFE.CT ct = fe.encrypt(setupParams.PK, revokeTree, m);
        LocalDateTime endEncrypt = LocalDateTime.now();

        // decrypt
        LocalDateTime startDecrypt = LocalDateTime.now();
        TraceRevokeFE.DecryptResult decryptResult = fe.decrypt(params, ct, SKX, id);
        LocalDateTime endDecrypt = LocalDateTime.now();

        // assert if correct
        Element gm = params.g.powZn(m);
        Element vrx = params.pairing.getZr().newZeroElement();
        List<Element> vr = fe.getVR(params, revokeTree);
        TraceRevokeFE.Theta theta = params.idTree.gamma.get(id);
        for (int i = 0; i < vr.size(); i++) {
            vrx = vrx.add(vr.get(i).mul(theta.theta.get(i)));
        }
        Element gmvrx = gm.powZn(vrx);
        assert decryptResult.decryptResult.equals(gmvrx);

        // time
        Duration setupDuration = Duration.between(startSetup, endSetup);
        Duration extractDuration = Duration.between(startKeyDer, endKeyDer);
        Duration encryptDuration = Duration.between(startEncrypt, endEncrypt);
        Duration decryptDuration = Duration.between(startDecrypt, endDecrypt);

        System.out.println("setupTime: " + setupDuration.toMillis());
        System.out.println("extractTime: " + extractDuration.toMillis());
        System.out.println("encryptTime: " + encryptDuration.toMillis());
        System.out.println("decryptTime: " + decryptDuration.toMillis());
    }

}