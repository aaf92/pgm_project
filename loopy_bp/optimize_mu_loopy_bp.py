import numpy as np
import pandas as pd
import networkx as nx
import sys
import pickle


# get prior network
def get_prior_network(net_file):
    # read in network
    prior_network = pd.read_csv(net_file, sep="\t")
    df = prior_network[["Gene_A", "Gene_B"]]
    df_cleaned = df.dropna()

    # Initialize graph (using NetworkX for graph structure)
    G = nx.from_pandas_edgelist(df_cleaned, source="Gene_A", target="Gene_B")
    return G

# factor potentital between nodes: this include lambda and mu parameter to add weight to prior edge and node state compatibility
def factor_potential(X_i, X_j, prior_edge, lambda_param, mu_param):
    if prior_edge == 1:
        # Higher potential if states are compatible (both active or both inactive)
        return np.exp(lambda_param + mu_param * (X_i == X_j))
    else:
        # Lower potential as no edge exists in the prior network
        return 1.0

# update message function
def update_message(i, j, G, messages, lambda_param, mu_param, epsilon=1e-6):
    # Calculate message M_{i -> j}(X_j) for X_j = 0 and X_j = 1
    new_message = {}
    for X_j in [0, 1]:
        # Sum over possible states of X_i
        message_sum = 0
        for X_i in [0, 1]:
            # Product of incoming messages from neighbors of i, excluding j
            incoming_product = 1
            for k in G.neighbors(i):
                if k != j:
                    incoming_product *= messages[(k, i)][X_i]
            # Calculate factor potential for the current X_i and X_j
            prior_edge = G.has_edge(i, j)
            potential = factor_potential(X_i, X_j, prior_edge, lambda_param, mu_param)
            message_sum += potential * incoming_product
        # Update message value for M_{i -> j}(X_j)
        new_message[X_j] = message_sum + epsilon

    # Normalize the message
    total = new_message[0] + new_message[1]
    new_message[0] /= total
    new_message[1] /= total
    return new_message


# initialize beliefs based on data derived activation state
def initialize_beliefs(G, observed_states):
    beliefs = {}
    for i in G.nodes():
        if i in observed_states:
            # If the node state is observed, set belief to the observed state
            state = abs(observed_states[i])
            beliefs[i] = {0: 0.0, 1: 0.0}
            beliefs[i][state] = 1.0
        else:
            # For unobserved nodes, initialize with equal probabilities
            beliefs[i] = {0: 0.5, 1: 0.5}
    return beliefs

# initialize messages based on data derived activation state
def initialize_messages(G, observed_states):
    messages = {}
    for (i, j) in G.edges():
        if i in observed_states:
            # If node i has an observed state, initialize message to reflect that state
            observed_state = abs(observed_states[i])
            messages[(i, j)] = {0: 0.0, 1: 0.0}
            messages[(i, j)][observed_state] = 1.0  # Strong message in favor of the observed state
        else:
            # For unobserved nodes, initialize messages to equal probabilities
            messages[(i, j)] = {0: 0.5, 1: 0.5}
            
        if j in observed_states:
            # Similarly, if node j has an observed state, initialize message to reflect that state
            observed_state = observed_states[j]
            messages[(j, i)] = {0: 0.0, 1: 0.0}
            messages[(j, i)][observed_state] = 1.0
        else:
            # For unobserved nodes, initialize messages to equal probabilities
            messages[(j, i)] = {0: 0.5, 1: 0.5}
    return messages


def loopy_bp(G, observed_states, T, lambda_param, mu_param):
    # Initialize beliefs and messages based on observed states
    messages = initialize_messages(G, observed_states)

    # Main loopy belief propagation loop
    for iteration in range(T):
        new_messages = {}
        for (i, j) in G.edges():
            new_message_ij = update_message(i, j, G, messages, lambda_param, mu_param)
            new_message_ji = update_message(j, i, G, messages, lambda_param, mu_param)

            new_messages[(i, j)] = new_message_ij
            new_messages[(j, i)] = new_message_ji
        messages.update(new_messages)

    # Calculate marginal beliefs using messages
    beliefs = {}
    for i in G.nodes():
        belief_i = {0: 1, 1: 1}
        for X_i in [0, 1]:
            for k in G.neighbors(i):
                belief_i[X_i] *= messages[(k, i)][X_i]
        # Normalize beliefs
        total = belief_i[0] + belief_i[1]
        if total == 0:
            beliefs[i] = {0: 0.5, 1: 0.5}
        else:
            beliefs[i] = {0: belief_i[0] / total, 1: belief_i[1] / total}
    return beliefs


# get posterior edge_df
def get_edge_df(beliefs, G):
    # estimate for posterior edges
    posterior_edges = {}
    for (i, j) in G.edges():
        # Compatibility-based posterior probability
        posterior_edges[(i, j)] = beliefs[i][1] * beliefs[j][1] + beliefs[i][0] * beliefs[j][0]

    edge_data = [{'Gene_A': i, 'Gene_B': j, 'posterior_edge': prob} 
                for (i, j), prob in posterior_edges.items()]

    # Create DataFrame from the list
    posterior_df = pd.DataFrame(edge_data, columns=['Gene_A', 'Gene_B', 'posterior_edge'])
    return posterior_df


def test_mu(posterior_G, moi_dict):
    overlap_score = 0
    for goi, moi in moi_dict.items():
        if len(moi) == 0:
            continue  # Skip empty MOI lists

        if goi in posterior_G:
            first_neighbors = set(posterior_G.neighbors(goi))
            intersect_moi = first_neighbors & set(moi)
            overlap_score += len(intersect_moi) / len(moi)
    return overlap_score


if __name__=="__main__":
    # inputs:
    net_file = sys.argv[1]
    states_file = sys.argv[2]
    output = sys.argv[3]
    prob_thresh = 0.70
    G = get_prior_network(net_file)

    observed_states = pd.read_csv(states_file)
    observed_states = observed_states.set_index('geneSymbol')['final_state'].to_dict()

    with open("/ix/djishnu/Aaron_F/PGM_project/loopy_bp/moi_dict.pkl", "rb") as file:
        moi_dict = pickle.load(file)


    # parameters:
    lambda_param = 0.6
    mu_params = np.linspace(0.05, 2, 20)
    T = 10000

    # test each mu_param for control module overlap scoring
    mu_score = {}
    for mu_param in mu_params:
        # Step1: get beliefs using current mu_param
        beliefs = loopy_bp(G, observed_states, T, lambda_param, mu_param)
        # Step2: transform beliefs into edge probabilities, then subset based on prob_threshold
        edge_data = get_edge_df(beliefs, G)
        edge_data = edge_data.loc[edge_data["posterior_edge"]>=prob_thresh,]
        cur_G = nx.from_pandas_edgelist(edge_data, source="Gene_A", target="Gene_B")
        # Step3: check for control module overlap and assign score to current mu_param
        overlap_score = test_mu(cur_G, moi_dict)
        mu_score[mu_param] = overlap_score
    
    # save
    score_df = pd.DataFrame(list(mu_score.items()), columns=["mu", "score"])
    score_df.to_csv(f"{output}/20241204_mu_score.csv", index=False)

