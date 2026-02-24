
# Correntes de Foucault — Simulação e Comparação de Materiais

Este repositório contém uma **simulação numérica de correntes de Foucault (eddy currents)** induzidas em diferentes materiais quando submetidos a um **campo magnético alternado**. O projeto compara alumínio, cobre e ferro em termos de:

- Penetração do campo magnético (skin effect)
- Densidade de corrente induzida \|J\|
- Potência dissipada por efeito Joule
- Defasagem temporal e atenuação do campo interno

Os dados gerados são utilizados para produzir gráficos comparativos e perfis de atenuação, além de uma tabela consolidada de parâmetros físicos.

---

## 📁 Estrutura do projeto

```

foucault/
├── main.py                # Script principal de simulação
├── data/                  # (Opcional) Dados gerados ou de entrada
├── outputs/               # Gráficos gerados
├── notebooks/             # Jupyter notebooks de análise
└── README.md              # (Este arquivo)

````

---

## 🧪 Requisitos

Este projeto foi desenvolvido em **Python** e depende de bibliotecas científicas como:

- `numpy`
- `scipy`
- `matplotlib`
- `pandas`

Instale dependências com:

```bash
pip install -r requirements.txt
````

> Se não existir `requirements.txt`, basta instalar as bibliotecas acima separadamente.

---

## 🚀 Como usar

1. Clone o repositório:

```bash
git clone https://github.com/vitor-souza-ime/foucault.git
cd foucault
```

2. Execute a simulação principal:

```bash
python main.py
```

3. Os gráficos e resultados serão gravados na pasta `outputs/` (se configurada no script).

---

## 📊 O que o código faz

A simulação:

* Define propriedades elétricas e magnéticas de materiais (σ, μr)
* Aplica um campo magnético alternado de 60 Hz
* Resolve a equação de **difusão magnética**
* Calcula:

  * Profundidade de penetração (*skin depth*)
  * Distribuição de correntes de Foucault
  * Potência dissipada
  * Perfil de atenuação e defasagem temporal

Esses resultados são usados para gerar:

* Mapas de corrente |J| dentro do material
* Mapas de densidade de potência dissipada
* Gráficos comparativos da atenuação do campo interno
* Curvas temporais do campo em profundidade fixa

---

## 📈 Exemplos de resultados

São gerados gráficos como:

* **Distribuição de correntes induzidas para cada material**
* **Potência dissipada no volume**
* **Perfil de atenuação (B(x))** — Simulado vs. Teórico
* **Evolução temporal do campo B em um ponto fixo**

A simulação produz também uma tabela comparativa de grandezas como:

| Grandeza              | Alumínio | Cobre  | Ferro  |
| --------------------- | -------- | ------ | ------ |
| Condutividade σ [S/m] | 3.50e7   | 5.80e7 | 1.00e7 |
| Permeabilidade μr     | 1        | 1      | 1000   |
| Skin depth δ [mm]     | ~10.98   | ~8.53  | ~0.65  |
| Defasagem             | …        | …      | …      |
| |J| máximo            | …        | …      | …      |
| Potência dissipada    | …        | …      | …      |

Estes valores mostram como o efeito de *skin* se torna mais forte em materiais de maior permeabilidade, limitando a penetração do campo e concentrando as correntes na superfície.

---

## 🧠 Interpretação física (resumo)

* **Materiais com maior condutividade** apresentam correntes induzidas mais intensas.
* **Materiais com maior permeabilidade (ferro)** exibem forte atenuação do campo, com skin depth muito pequena.
* **Fenômeno de skin effect** é bem evidenciado nos gráficos, comparando B(x) simulado e teórico.
* A defasagem entre campo aplicado e resposta interna depende da difusão magnética do material.

---

## 📝 Referências

Esta simulação está baseada em princípios eletromagnéticos clássicos da indução, efeito pelicular e equação de difusão para campos magnéticos em meios condutores:

* Lei de Faraday
* Lei de Ohm
* Equação de difusão magnética

---

## 👤 Autor

**Vítor Amadeu Souza**
Simulações e visualizações de correntes de Foucault.


