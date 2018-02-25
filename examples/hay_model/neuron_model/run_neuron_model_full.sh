SEED=1

while [  $SEED -lt 100 ]; do
  python simulate_hay_model_store_current_from_all_sections.py $SEED
  python combine_seed_data.py $SEED
  SEED=$((SEED+10))
done

python combine_different_seeds_to_single_neuron_data.py 1
