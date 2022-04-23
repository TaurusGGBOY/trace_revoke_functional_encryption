package algorithm;

import it.unisa.dia.gas.jpbc.*;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;
import it.unisa.dia.gas.plaf.jpbc.pairing.a.TypeACurveGenerator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.*;

public class TraceRevokeFE {

    public static final String paramFilePath = "params.properties";
    public static final String propertyFilePath = "properties.properties";

    class Params {
        Pairing pairing;
        Field G;
        Element g;
        Element h;
        BigInteger L;
        IdTree idTree;

        public Params(Pairing pairing, Field g, Element g1, Element h, BigInteger l, IdTree idTree) {
            this.pairing = pairing;
            G = g;
            this.g = g1;
            this.h = h;
            L = l;
            this.idTree = idTree;
        }
    }

    class PK {
        Params params;
        List<Element> hs;

        public PK(Params params, List<Element> hs) {
            this.params = params;
            this.hs = hs;
        }
    }

    class MSK {

        List<Element> ss;
        List<Element> ts;

        public MSK(List<Element> ss, List<Element> ts) {
            this.ss = ss;
            this.ts = ts;
        }
    }

    class SKX {
        Element xs;
        Element xt;

        public SKX(Element xs, Element xt) {
            this.xs = xs;
            this.xt = xt;
        }
    }

    class IdTree {
        List<Theta> gamma;

        public IdTree(List<Theta> gamma) {
            this.gamma = gamma;
        }
    }

    class Theta {
        List<Element> theta;

        public Theta(List<Element> theta) {
            this.theta = theta;
        }
    }

    class SetupParams {
        public SetupParams(PK PK, MSK msk) {
            this.PK = PK;
            this.msk = msk;
        }

        PK PK;
        MSK msk;
    }

    class CT {
        Element d0;
        Element d1;
        List<Element> es;
        IdTree revokeTree;

        public CT(Element d0, Element d1, List<Element> es, IdTree revokeTree) {
            this.d0 = d0;
            this.d1 = d1;
            this.es = es;
            this.revokeTree = revokeTree;
        }
    }

    class DecryptResult {
        public DecryptResult(Element decryptResult) {
            this.decryptResult = decryptResult;
        }

        Element decryptResult;
    }

    public void curveInit(int lambda) throws IOException {
        // curve 256bit, limit field 512bit
        // A E A1 F four curves
        PairingParametersGenerator generator = new TypeACurveGenerator(lambda, lambda * 2);
        PairingParameters parameters = generator.generate();
        BufferedWriter writer = new BufferedWriter(new FileWriter(paramFilePath));
        writer.write(parameters.toString());
        writer.close();
    }

    public List<Element> getVR(Params params, IdTree revokeTree){
        Gauss gauss = new Gauss(revokeTree, params.pairing);
        return gauss.getOneSolutionVector();
    }

    public SetupParams setup(int lambda, int t, int r) throws IOException {
        // init curve params
        LocalDateTime curvInitStart = LocalDateTime.now();
        File paramFile = new File(paramFilePath);
        if (!paramFile.exists() || paramFile.length() == 0) {
            curveInit(lambda);
        }



        // get params
        Pairing pairing = PairingFactory.getPairing(paramFilePath);
        LocalDateTime curvInitEnd = LocalDateTime.now();
        Duration curvInitDuration = Duration.between(curvInitStart, curvInitEnd);
        System.out.println("curvInit time: " + curvInitDuration.toMillis());

        Field G = pairing.getG1();
        BigInteger L = G.getOrder();
        Element g = G.newRandomElement().getImmutable();
        Element h = G.newRandomElement().getImmutable();
        int ell = t + r + 1;

        // init s t and h_i
        List<Element> ss = new ArrayList<>();
        List<Element> ts = new ArrayList<>();
        List<Element> hs = new ArrayList<>();
        for (int i = 0; i < ell; i++) {
            ss.add(pairing.getZr().newRandomElement().getImmutable());
            ts.add(pairing.getZr().newRandomElement().getImmutable());
            Element gs = g.powZn(ss.get(i)).getImmutable();
            Element ht = h.powZn(ts.get(i)).getImmutable();
            hs.add(gs.mul(ht).getImmutable());
        }

        // init gamma
        List<TraceRevokeFE.Theta> thetas = new ArrayList<>();
        // user number must less than ell so that VR can be solved
        for (int i = 0; i < ell; i++) {
            List<Element> theta = new ArrayList<>();
            for (int j = 0; j < ell; j++) {
                theta.add(pairing.getZr().newRandomElement().getImmutable());
            }
            thetas.add(new Theta(theta));
        }
        IdTree idTree = new IdTree(thetas);

        // init params
        Params params = new Params(pairing, G, g, h, L, idTree);

        // init mpk and msk
        PK pk = new PK(params, hs);
        MSK msk = new MSK(ss, ts);

        // init setup params
        SetupParams setupParams = new SetupParams(pk, msk);

        return setupParams;
    }

    public SKX keyGen(Params params, int id, MSK msk) {
        // get thetaID
        Theta thetaID = params.idTree.gamma.get(id);

        // get s and t
        List<Element> ss = msk.ss;
        List<Element> ts = msk.ts;

        // compute intermedia result
        Element xs = params.pairing.getZr().newZeroElement();
        Element xt = params.pairing.getZr().newZeroElement();
        for (int i = 0; i < thetaID.theta.size(); i++) {
            xs = xs.add(ss.get(i).mul(thetaID.theta.get(i)));
            xt = xt.add(ts.get(i).mul(thetaID.theta.get(i)));
        }

        return new SKX(xs.getImmutable(), xt.getImmutable());
    }

    public CT encrypt(PK pk, IdTree revokeTree,Element m) {
        Params params = pk.params;

        // compute vr
        List<Element> vr = getVR(params, revokeTree);

        // random choose
        Element r = params.pairing.getZr().newRandomElement().getImmutable();

        // compute d0 d1 ei
        // compute y
        Element d0 = params.g.powZn(r).getImmutable();
        Element d1 = params.h.powZn(r).getImmutable();
        List<Element> ys = new ArrayList<>();
        List<Element> es = new ArrayList<>();
        for (int i = 0; i < vr.size(); i++) {
            ys.add(vr.get(i).mul(m).getImmutable());
            Element gy = params.g.powZn(ys.get(i));
            Element hr = pk.hs.get(i).powZn(r);
            es.add(gy.mul(hr).getImmutable());
        }

        return new CT(d0, d1, es, revokeTree);
    }

    public DecryptResult decrypt(Params params, CT ct, SKX skx, int id) {
        // get thetaid
        Theta thetaID = params.idTree.gamma.get(id);

        // compute upper part
        Element decryptResult = params.pairing.getG1().newOneElement();
        for (int i = 0; i < thetaID.theta.size(); i++) {
            Element ex = ct.es.get(i).powZn(thetaID.theta.get(i));
            decryptResult = decryptResult.mul(ex);
        }

        // div lower part
        Element ds = ct.d0.powZn(skx.xs).getImmutable();
        Element dt = ct.d1.powZn(skx.xt).getImmutable();
        decryptResult = decryptResult.div(ds.mul(dt));
        return new DecryptResult(decryptResult);
    }
}