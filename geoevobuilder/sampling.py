import datetime
import os
import random

import numpy as np

# from geoevobuilder.distributionR.samplingR import samplingr
import torch
import torch.nn.functional as F

from geoevobuilder import (
    AA,
    AA_reverse,
    Val_builder,
    batch_fea,
    batch_tri,
    fetchseq,
)

standard_aa_names = {
    "ALA": 0,
    "CYS": 1,
    "ASP": 2,
    "GLU": 3,
    "PHE": 4,
    "GLY": 5,
    "HIS": 6,
    "ILE": 7,
    "LYS": 8,
    "LEU": 9,
    "MET": 10,
    "ASN": 11,
    "PRO": 12,
    "GLN": 13,
    "ARG": 14,
    "SER": 15,
    "THR": 16,
    "VAL": 17,
    "TRP": 18,
    "TYR": 19,
}
standard_aa_names_r = dict([val, key] for key, val in standard_aa_names.items())


def samPling(args, models):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Val_builder(
        os.path.join(args.inputPATH, args.inputfile), args.chainID, args.inputfile[:-4] + "_Val.pdb"
    )
    seq_alpha, seqnumber = fetchseq(os.path.join(args.inputPATH, args.inputfile), args.chainID)
    protein_graph = batch_fea(args.inputfile[:-4] + "_Val.pdb", args.chainID, len(seq_alpha))
    # features
    y = torch.tensor(seqnumber).to(device)

    x = torch.from_numpy(protein_graph["Node_features"]).float().to(device)
    edge_index = torch.from_numpy(protein_graph["edge_index"]).T.long().to(device)
    edge_attr = torch.from_numpy(protein_graph["edge_attrs"]).float().to(device)

    triang = batch_tri(
        args.inputfile[:-4] + "_Val.pdb", args.chainID, len(seq_alpha), protein_graph["edge_index"]
    )

    if args.purpose == 0:
        model_se = [models["model%d" % i].to(device).eval() for i in [1, 2]]
    elif args.purpose == 1:
        model_se = models["model_MSA"].to(device).eval()
        MSA_info = torch.load(args.MSA).to(device)
    os.system(f"rm {args.inputfile[:-4]}_Val.pdb")
    print("\n" + "-" * 50)
    with torch.no_grad():
        average_acc = 0
        final_seqs = open(
            os.path.join(args.outputPATH, args.inputfile[:-4])
            + "_%smasked_" % (str(args.SM))
            + str(datetime.date.today())
            + ".fa",
            "w",
        )
        for output_iter in range(args.n):
            num = range(20)
            np.random.seed(output_iter)
            if args.SM == 100:
                y_initial = (
                    torch
                    .from_numpy(np.random.choice(num, size=x.size(0), replace=True))
                    .long()
                    .to(device)
                )
            elif args.SM == 0:
                y_initial = y
            else:
                y_ = y.clone()
                mask_pos = np.random.choice(
                    [i for i in range(x.size(0))],
                    size=int(x.size(0) * args.SM / 100),
                    replace=False,
                )
                y_initial = y_
                # y_initial[mask_pos] = 0 ##residues at all masked positions are replaced using Ala
                y_initial[mask_pos] = (
                    torch
                    .from_numpy(
                        np.random.choice(
                            [i for i in range(20)],
                            size=int(x.size(0) * args.SM / 100),
                            replace=True,
                        )
                    )
                    .unsqueeze(0)
                    .to(device)
                )  # random mutations for these masked positions
            condition = torch.scatter(
                torch.zeros(y_initial.size(0), 20).to(device), 1, y_initial.unsqueeze(1), 1
            ).float()
            seq_initial = "".join(AA_reverse[index.item()] for index in y_initial)
            batch_tokens, decode_tockens = esm2emb(seq_initial)
            esm2in = torch.nn.functional.one_hot(decode_tockens.detach(), num_classes=33).to(device)

            if args.Fixed:
                fix = args.Fixed  # .split(',')
                fixid = [i.split("_") for i in fix]
                fixid_list = [int(i) for (i, j, k) in fixid]
                for i, j, k in fixid:
                    if j != args.chainID:
                        raise ValueError(
                            "--Fixed, %s dose not mathed the designed chain %s" % (j, args.chainID)
                        )
                    resid, fixres = int(i), k
                    condition[resid][:20] = F.one_hot(torch.tensor(AA[fixres], device=device), 20)
            aa = ""
            proba = 0
            sim_old_new = []
            seq = []
            if args.ST > 0.5:
                ks = [3] * int(condition.size(0) * 1.5) + [1] * int(condition.size(0) / 2)
            else:
                ks = [3] * int(condition.size(0) / 2) + [1] * int(condition.size(0) / 2)

            for i, k in enumerate(ks):
                if args.purpose == 0:
                    pred1 = model_se[0](
                        x=x,
                        edge_index=edge_index,
                        edge_attr=edge_attr.float(),
                        condition=condition,
                        triangle_nodes=torch
                        .from_numpy(triang["triangleindex"])
                        .T.long()
                        .to(device),
                        triangle_edges=torch
                        .from_numpy(triang["triangle2line"])
                        .T.long()
                        .to(device),
                        esm2in=esm2in.float(),
                    )
                    pred2 = model_se[0](
                        x=x,
                        edge_index=edge_index,
                        edge_attr=edge_attr.float(),
                        condition=condition,
                        triangle_nodes=torch
                        .from_numpy(triang["triangleindex"])
                        .T.long()
                        .to(device),
                        triangle_edges=torch
                        .from_numpy(triang["triangle2line"])
                        .T.long()
                        .to(device),
                        esm2in=esm2in.float(),
                    )
                    pred1 = (pred1 + pred2) / 2
                    proba_softmax = F.softmax(pred1 / args.ST, 1)
                else:
                    pred1 = model_se(
                        x=x,
                        edge_index=edge_index,
                        edge_attr=edge_attr.float(),
                        condition=condition,
                        triangle_nodes=torch
                        .from_numpy(triang["triangleindex"])
                        .T.long()
                        .to(device),
                        triangle_edges=torch
                        .from_numpy(triang["triangle2line"])
                        .T.long()
                        .to(device),
                        esm2in=esm2in.float(),
                        MSA=MSA_info,
                    )
                    proba_softmax = F.softmax(pred1 / args.ST, 1)
                non_optimal = (
                    (condition[:, :20].argmax(1) - proba_softmax.argmax(1))
                    .nonzero()
                    .squeeze(dim=1)
                    .tolist()
                )
                if args.Fixed:
                    non_optimal = [i for i in non_optimal if i not in fixid_list]
                muteres = random.sample(non_optimal, min(len(non_optimal), k))
                step_acc = sum(condition[:, :20].argmax(1) == pred1.argmax(1)).item() / len(y)
                sim_old_new.append(step_acc)
                seq.append(
                    "".join(AA_reverse[index.item()] for index in condition[:, :20].argmax(1))
                )
                if (
                    len(sim_old_new) > 10
                    and sim_old_new[-1] - sim_old_new[-10] < 1e-6
                    and step_acc == 1
                ) or len(non_optimal) == 0:
                    break
                for j in muteres:
                    condition[j] = F.one_hot(
                        torch.multinomial(proba_softmax[j], num_samples=1), 20
                    ).float()
                seqs_changes = "".join(AA_reverse[index.item()] for index in condition.argmax(1))
                batch_tokens, decode_tockens = esm2emb(seqs_changes)
                esm2in = torch.nn.functional.one_hot(decode_tockens.detach(), num_classes=33).to(
                    device
                )

                logP = torch.log(
                    torch.gather(
                        F.softmax(pred1, 1), index=condition[:, :20].argmax(1).unsqueeze(1), dim=1
                    )
                )
                logits = (
                    torch
                    .gather(pred1, dim=1, index=condition[:, :20].argmax(1).unsqueeze(1))
                    .mean()
                    .item()
                )
                if args.verbose == 1:
                    print(
                        seq[-1],
                        "convergent identity: %.4f," % step_acc,
                        "iterative steps: %d," % i,
                        "#(unconvergent residues): %d," % len(non_optimal),
                        "logP: %.4f," % logP.mean().item(),
                        "logits: %.4f" % logits,
                    )
            seq = "".join(AA_reverse[index.item()] for index in condition[:, :20].argmax(1))
            acc = sum(condition[:, :20].argmax(1) == y).item() / len(y)
            print(
                "Convergent sequence of Epoch %d," % output_iter
                + " " * (3 - len(str(output_iter))),
                seq + ",",
                "relative to native identity: %.4f," % acc,
                "logP: %.4f,"
                % (torch.log(F.softmax(pred1.detach(), 1).max(1)[0] + 1e-10).mean().item()),
                "logits: %.4f" % (pred1.detach().max(1)[0].mean().item()),
            )
            final_seqs.write(
                ">%d" % output_iter
                + " | convergent identity: %.4f" % step_acc
                + " | iterative steps: %d" % i
                + " | #(unconvergent residues): %d" % len(non_optimal)
                + " | relative to native identity: %.4f" % acc
                + " | logP: %.4f" % torch.log(F.softmax(pred1, 1).max(1)[0] + 1e-10).mean().item()
                + " | logits: %.4f" % pred1.max(1)[0].mean().item()
                + "\n"
            )
            final_seqs.write(seq + "\n")
            average_acc += acc
        final_seqs.close()
        print(f"Average sequence recovery: {average_acc / (output_iter + 1):.4f}")
        print(
            "Writing sequences to "
            + os.path.abspath(
                os.path.join(args.outputPATH, args.inputfile[:-4])
                + "_%smasked_" % (str(args.SM))
                + str(datetime.date.today())
                + ".fa"
            )
        )
        print("Done!")
        return


import esm

# Load ESM-2 model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
esm2model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
batch_converter = alphabet.get_batch_converter()
esm2model.to(device).eval()


def esm2emb(seq):
    data = [("tmp", seq)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = esm2model(batch_tokens.to(device), repr_layers=[36], return_contacts=True)
        return batch_tokens[0][1:-1], results["logits"].argmax(-1)[0][1:-1]
