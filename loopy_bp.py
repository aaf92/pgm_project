import numpy as np
import pandas as pd
import networkx as nx
import sys


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


if __name__=="__main__":
    # inputs:
    net_file = sys.argv[1]
    states_file = sys.argv[2]
    output = sys.argv[3]
    G = get_prior_network(net_file)

    observed_states = pd.read_csv(states_file)
    observed_states = observed_states.set_index('geneSymbol')['final_state'].to_dict()


    # parameters:
    lambda_param = 0.6
    mu_param = 0.8
    T = 1000

    # Run Loopy belief prop:
    # Initialize beliefs and messages based on observed states
    messages = initialize_messages(G, observed_states)
    beliefs = initialize_beliefs(G, observed_states)
    
    
    # main for loopy bp
    for iteration in range(T):
        new_messages = {}
        # Update messages for each edge
        for (i, j) in G.edges():
            new_message_ij = update_message(i, j, G, messages, lambda_param, mu_param)
            new_message_ji = update_message(j, i, G, messages, lambda_param, mu_param)
            
            # Store updated messages
            new_messages[(i, j)] = new_message_ij
            new_messages[(j, i)] = new_message_ji
        # Update all messages
        messages.update(new_messages)
    
    # now calculate marginal beliefs using messages
    beliefs = {}
    for i in G.nodes():
        belief_i = {0: 1, 1: 1}
        for X_i in [0, 1]:
            for k in G.neighbors(i):
                belief_i[X_i] *= messages[(k, i)][X_i]
        # Normalize the beliefs
        total = belief_i[0] + belief_i[1]
        beliefs[i] = {0: belief_i[0] / total, 1: belief_i[1] / total}

    # estimate for posterior edges
    posterior_edges = {}
    for (i, j) in G.edges():
        # Compatibility-based posterior probability
        posterior_edges[(i, j)] = beliefs[i][1] * beliefs[j][1] + beliefs[i][0] * beliefs[j][0]

    # Display posterior probabilities of edges
    for edge, posterior in posterior_edges.items():
        print(f"Edge {edge}: Posterior Probability = {posterior:.4f}")

    edge_data = [{'Gene_A': i, 'Gene_B': j, 'posterior_edge': prob} 
                 for (i, j), prob in posterior_edges.items()]

    # Create DataFrame from the list
    posterior_df = pd.DataFrame(edge_data, columns=['Gene_A', 'Gene_B', 'posterior_edge'])
    # save
    posterior_df.to_csv(f"{output}/edge_data.csv", index=False)
